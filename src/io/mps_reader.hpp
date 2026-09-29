#pragma once
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cctype>
#include <cstdio>
#include "../model/model.hpp"

namespace sov {

class MpsReader {
public:
    static Status read(const std::string& path, Model& out_model, std::string& error_msg) {
        std::ifstream in(path);
        if (!in.is_open()) {
            error_msg = "Cannot open MPS file: " + path;
            return Status::InputError;
        }

        Model m;
        m.name = "mps_model";
        std::string obj_row_name = "";
        bool has_obj_sense_Wait = false;

        enum Section {
            SEC_NONE, SEC_NAME, SEC_OBJSENSE, SEC_OBJNAME,
            SEC_ROWS, SEC_COLUMNS, SEC_RHS, SEC_RANGES,
            SEC_BOUNDS, SEC_QUADOBJ, SEC_QMATRIX, SEC_SOS, SEC_ENDATA
        } sec = SEC_NONE;

        std::unordered_map<std::string, Idx> row_map;
        std::unordered_map<std::string, Idx> col_map;
        std::vector<RowType> row_types;
        std::vector<double> row_rhs;
        std::vector<double> row_range;
        std::vector<bool> row_has_range;
        std::vector<Triplet> triplets;
        std::vector<Triplet> q_triplets;
        bool in_integer_block = false;

        auto ensure_col = [&](const std::string& cname) -> Idx {
            auto it = col_map.find(cname);
            if (it != col_map.end()) return it->second;
            Idx j = m.ncol++;
            col_map[cname] = j;
            m.col_names.push_back(cname);
            m.cost.push_back(0.0);
            m.collo.push_back(0.0);
            m.colup.push_back(kInf);
            m.vartype.push_back(in_integer_block ? VarType::Integer : VarType::Continuous);
            return j;
        };

        auto parse_double = [](std::string s) -> double {
            for (char& ch : s) {
                if (ch == 'D' || ch == 'd') ch = 'E';
            }
            return std::stod(s);
        };

        std::string line;
        int line_num = 0;
        while (std::getline(in, line)) {
            ++line_num;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty() || line[0] == '*') continue;

            // Trim trailing comment or whitespace
            std::istringstream iss(line);
            std::vector<std::string> tok;
            std::string t;
            while (iss >> t) tok.push_back(t);
            if (tok.empty()) continue;

            // Section header if column 0 is non-space or known keyword
            std::string first_up = tok[0];
            for (char& c : first_up) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));

            if (!std::isspace(static_cast<unsigned char>(line[0]))) {
                if (first_up == "NAME") {
                    if (tok.size() >= 2) m.name = tok[1];
                    sec = SEC_NAME;
                    continue;
                } else if (first_up == "OBJSENSE") {
                    if (tok.size() >= 2) {
                        std::string s = tok[1];
                        for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                        if (s.find("MAX") == 0) m.sense = Sense::Maximize;
                        else m.sense = Sense::Minimize;
                    } else {
                        has_obj_sense_Wait = true;
                        sec = SEC_OBJSENSE;
                    }
                    continue;
                } else if (first_up == "OBJNAME") {
                    if (tok.size() >= 2) obj_row_name = tok[1];
                    else sec = SEC_OBJNAME;
                    continue;
                } else if (first_up == "ROWS") { sec = SEC_ROWS; continue; }
                else if (first_up == "COLUMNS") { sec = SEC_COLUMNS; continue; }
                else if (first_up == "RHS") { sec = SEC_RHS; continue; }
                else if (first_up == "RANGES") { sec = SEC_RANGES; continue; }
                else if (first_up == "BOUNDS") { sec = SEC_BOUNDS; continue; }
                else if (first_up == "QUADOBJ" || first_up == "QSECTION") { sec = SEC_QUADOBJ; continue; }
                else if (first_up == "QMATRIX") { sec = SEC_QMATRIX; continue; }
                else if (first_up == "SOS") { sec = SEC_SOS; continue; }
                else if (first_up == "ENDATA") { sec = SEC_ENDATA; break; }
            }

            if (sec == SEC_OBJSENSE && has_obj_sense_Wait) {
                if (first_up.find("MAX") == 0) m.sense = Sense::Maximize;
                else m.sense = Sense::Minimize;
                has_obj_sense_Wait = false;
                continue;
            }
            if (sec == SEC_OBJNAME && obj_row_name.empty()) {
                obj_row_name = tok[0];
                continue;
            }

            if (sec == SEC_ROWS) {
                if (tok.size() < 2) continue;
                char rt = static_cast<char>(std::toupper(static_cast<unsigned char>(tok[0][0])));
                const std::string& rname = tok[1];
                if (rt == 'N') {
                    if (obj_row_name.empty()) {
                        obj_row_name = rname;
                    }
                    // Free rows beyond the objective row are ignored per 3.2.3
                } else {
                    if (row_map.count(rname)) {
                        error_msg = "Duplicate row name at line " + std::to_string(line_num);
                        return Status::InputError;
                    }
                    Idx i = m.nrow++;
                    row_map[rname] = i;
                    m.row_names.push_back(rname);
                    RowType type = RowType::E;
                    if (rt == 'L') type = RowType::L;
                    else if (rt == 'G') type = RowType::G;
                    else if (rt == 'E') type = RowType::E;
                    row_types.push_back(type);
                    row_rhs.push_back(0.0);
                    row_range.push_back(0.0);
                    row_has_range.push_back(false);
                }
            } else if (sec == SEC_COLUMNS) {
                // Check for MARKER
                bool is_marker = false;
                for (const auto& tk : tok) {
                    if (tk.find("INTORG") != std::string::npos) {
                        in_integer_block = true;
                        is_marker = true;
                    } else if (tk.find("INTEND") != std::string::npos) {
                        in_integer_block = false;
                        is_marker = true;
                    }
                }
                if (is_marker) continue;
                if (tok.size() < 3) continue;

                const std::string& cname = tok[0];
                Idx j = ensure_col(cname);
                for (size_t pos = 1; pos + 1 < tok.size(); pos += 2) {
                    const std::string& rname = tok[pos];
                    double val = parse_double(tok[pos + 1]);
                    if (rname == obj_row_name) {
                        m.cost[j] += val;
                    } else {
                        auto it = row_map.find(rname);
                        if (it != row_map.end()) {
                            triplets.push_back({it->second, j, val});
                        }
                    }
                }
            } else if (sec == SEC_RHS) {
                // In free format RHS may have or omit the RHS set name:
                // If tok.size() is odd (3 or 5), tok[0] is RHS set name; if even (2 or 4), tok[0] is row name
                size_t start_pos = (tok.size() % 2 == 1) ? 1 : 0;
                for (size_t pos = start_pos; pos + 1 < tok.size(); pos += 2) {
                    const std::string& rname = tok[pos];
                    double val = parse_double(tok[pos + 1]);
                    if (rname == obj_row_name) {
                        // Standard MPS convention: objective offset = -RHS
                        m.obj_offset = -val;
                    } else {
                        auto it = row_map.find(rname);
                        if (it != row_map.end()) {
                            row_rhs[it->second] = val;
                        }
                    }
                }
            } else if (sec == SEC_RANGES) {
                size_t start_pos = (tok.size() % 2 == 1) ? 1 : 0;
                for (size_t pos = start_pos; pos + 1 < tok.size(); pos += 2) {
                    const std::string& rname = tok[pos];
                    double val = parse_double(tok[pos + 1]);
                    auto it = row_map.find(rname);
                    if (it != row_map.end()) {
                        row_range[it->second] = val;
                        row_has_range[it->second] = true;
                    }
                }
            } else if (sec == SEC_BOUNDS) {
                if (tok.size() < 3) continue;
                std::string btype = tok[0];
                for (char& c : btype) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                // Format: <type> <bnd_name> <col_name> [<val>]  OR  <type> <col_name> [<val>]
                std::string cname;
                double val = 0.0;
                bool has_val = (btype == "LO" || btype == "UP" || btype == "FX" || btype == "LI" || btype == "UI" || btype == "SC");
                if (has_val) {
                    if (tok.size() >= 4) {
                        cname = tok[2];
                        val = parse_double(tok[3]);
                    } else {
                        cname = tok[1];
                        val = parse_double(tok[2]);
                    }
                } else {
                    cname = (tok.size() >= 3) ? tok[2] : tok[1];
                }
                auto it = col_map.find(cname);
                if (it == col_map.end()) {
                    error_msg = "Unknown column in BOUNDS: " + cname;
                    return Status::InputError;
                }
                Idx j = it->second;
                if (btype == "LO") {
                    m.collo[j] = val;
                } else if (btype == "UP") {
                    if (val < 0.0 && m.collo[j] == 0.0) {
                        m.collo[j] = -kInf; // Legacy negative UP rule (3.2.7)
                    }
                    m.colup[j] = val;
                } else if (btype == "FX") {
                    m.collo[j] = val;
                    m.colup[j] = val;
                } else if (btype == "FR") {
                    m.collo[j] = -kInf;
                    m.colup[j] = kInf;
                } else if (btype == "MI") {
                    m.collo[j] = -kInf;
                } else if (btype == "PL") {
                    m.colup[j] = kInf;
                } else if (btype == "BV") {
                    m.vartype[j] = VarType::Binary;
                    m.collo[j] = 0.0;
                    m.colup[j] = 1.0;
                } else if (btype == "LI") {
                    m.vartype[j] = VarType::Integer;
                    m.collo[j] = val;
                } else if (btype == "UI") {
                    m.vartype[j] = VarType::Integer;
                    m.colup[j] = val;
                }
            } else if (sec == SEC_QUADOBJ || sec == SEC_QMATRIX) {
                if (tok.size() < 3) continue;
                auto it1 = col_map.find(tok[0]);
                auto it2 = col_map.find(tok[1]);
                if (it1 == col_map.end() || it2 == col_map.end()) continue;
                Idx i = it1->second;
                Idx j = it2->second;
                double val = parse_double(tok[2]);
                if (sec == SEC_QMATRIX) {
                    if (i <= j) q_triplets.push_back({i, j, val});
                } else {
                    if (i > j) std::swap(i, j);
                    q_triplets.push_back({i, j, val});
                }
            }
        }

        // Assemble row bounds (with RANGES per 3.2.6)
        m.rowlo.resize(m.nrow);
        m.rowup.resize(m.nrow);
        for (Idx i = 0; i < m.nrow; ++i) {
            double b = row_rhs[i];
            if (!row_has_range[i]) {
                if (row_types[i] == RowType::L) { m.rowlo[i] = -kInf; m.rowup[i] = b; }
                else if (row_types[i] == RowType::G) { m.rowlo[i] = b; m.rowup[i] = kInf; }
                else { m.rowlo[i] = b; m.rowup[i] = b; }
            } else {
                double R = row_range[i];
                if (row_types[i] == RowType::G) {
                    m.rowlo[i] = b;
                    m.rowup[i] = b + std::abs(R);
                } else if (row_types[i] == RowType::L) {
                    m.rowlo[i] = b - std::abs(R);
                    m.rowup[i] = b;
                } else {
                    if (R >= 0.0) { m.rowlo[i] = b; m.rowup[i] = b + R; }
                    else { m.rowlo[i] = b + R; m.rowup[i] = b; }
                }
            }
        }

        m.A = CscMatrix::from_triplets(m.nrow, m.ncol, triplets);
        m.Q = CscMatrix::from_triplets(m.ncol, m.ncol, q_triplets);
        out_model = std::move(m);
        return out_model.validate(error_msg);
    }

    static bool write(const std::string& path, const Model& m) {
        std::ofstream out(path);
        if (!out.is_open()) return false;
        char buf[256];
        out << "NAME          " << m.name << "\n";
        out << "OBJSENSE\n  " << (m.sense == Sense::Maximize ? "MAX" : "MIN") << "\n";
        out << "ROWS\n";
        out << " N  OBJ\n";
        for (Idx i = 0; i < m.nrow; ++i) {
            std::string rname = (i < static_cast<Idx>(m.row_names.size())) ? m.row_names[i] : ("R" + std::to_string(i));
            char rtype = 'E';
            if (!is_finite_bound(m.rowlo[i]) && is_finite_bound(m.rowup[i])) rtype = 'L';
            else if (is_finite_bound(m.rowlo[i]) && !is_finite_bound(m.rowup[i])) rtype = 'G';
            else if (m.rowlo[i] < m.rowup[i]) rtype = 'G'; // Ranged handled via RANGES
            out << " " << rtype << "  " << rname << "\n";
        }
        out << "COLUMNS\n";
        bool in_int = false;
        for (Idx j = 0; j < m.ncol; ++j) {
            bool is_int = (m.vartype[j] == VarType::Integer || m.vartype[j] == VarType::Binary);
            if (is_int && !in_int) {
                out << "    MARK0000  'MARKER'                 'INTORG'\n";
                in_int = true;
            } else if (!is_int && in_int) {
                out << "    MARK0001  'MARKER'                 'INTEND'\n";
                in_int = false;
            }
            std::string cname = (j < static_cast<Idx>(m.col_names.size())) ? m.col_names[j] : ("C" + std::to_string(j));
            if (m.cost[j] != 0.0) {
                std::snprintf(buf, sizeof(buf), "    %-10s %-10s %.17g\n", cname.c_str(), "OBJ", m.cost[j]);
                out << buf;
            }
            for (Off k = m.A.start[j]; k < m.A.start[j + 1]; ++k) {
                Idx i = m.A.index[k];
                std::string rname = (i < static_cast<Idx>(m.row_names.size())) ? m.row_names[i] : ("R" + std::to_string(i));
                std::snprintf(buf, sizeof(buf), "    %-10s %-10s %.17g\n", cname.c_str(), rname.c_str(), m.A.value[k]);
                out << buf;
            }
        }
        if (in_int) {
            out << "    MARK0001  'MARKER'                 'INTEND'\n";
        }
        out << "RHS\n";
        if (m.obj_offset != 0.0) {
            std::snprintf(buf, sizeof(buf), "    RHS1       %-10s %.17g\n", "OBJ", -m.obj_offset);
            out << buf;
        }
        for (Idx i = 0; i < m.nrow; ++i) {
            std::string rname = (i < static_cast<Idx>(m.row_names.size())) ? m.row_names[i] : ("R" + std::to_string(i));
            double rhs = is_finite_bound(m.rowlo[i]) ? m.rowlo[i] : m.rowup[i];
            if (rhs != 0.0 && is_finite_bound(rhs)) {
                std::snprintf(buf, sizeof(buf), "    RHS1       %-10s %.17g\n", rname.c_str(), rhs);
                out << buf;
            }
        }
        out << "RANGES\n";
        for (Idx i = 0; i < m.nrow; ++i) {
            if (is_finite_bound(m.rowlo[i]) && is_finite_bound(m.rowup[i]) && m.rowlo[i] < m.rowup[i]) {
                std::string rname = (i < static_cast<Idx>(m.row_names.size())) ? m.row_names[i] : ("R" + std::to_string(i));
                std::snprintf(buf, sizeof(buf), "    RNG1       %-10s %.17g\n", rname.c_str(), m.rowup[i] - m.rowlo[i]);
                out << buf;
            }
        }
        out << "BOUNDS\n";
        for (Idx j = 0; j < m.ncol; ++j) {
            std::string cname = (j < static_cast<Idx>(m.col_names.size())) ? m.col_names[j] : ("C" + std::to_string(j));
            if (m.vartype[j] == VarType::Binary && m.collo[j] == 0.0 && m.colup[j] == 1.0) {
                std::snprintf(buf, sizeof(buf), " BV BND1       %s\n", cname.c_str());
                out << buf;
                continue;
            }
            if (!is_finite_bound(m.collo[j]) && !is_finite_bound(m.colup[j])) {
                std::snprintf(buf, sizeof(buf), " FR BND1       %s\n", cname.c_str());
                out << buf;
                continue;
            }
            if (is_finite_bound(m.collo[j]) && is_finite_bound(m.colup[j]) && m.collo[j] == m.colup[j]) {
                std::snprintf(buf, sizeof(buf), " FX BND1       %-10s %.17g\n", cname.c_str(), m.collo[j]);
                out << buf;
                continue;
            }
            if (!is_finite_bound(m.collo[j])) {
                std::snprintf(buf, sizeof(buf), " MI BND1       %s\n", cname.c_str());
                out << buf;
            } else if (m.collo[j] != 0.0) {
                std::snprintf(buf, sizeof(buf), " LO BND1       %-10s %.17g\n", cname.c_str(), m.collo[j]);
                out << buf;
            }
            if (is_finite_bound(m.colup[j])) {
                std::snprintf(buf, sizeof(buf), " UP BND1       %-10s %.17g\n", cname.c_str(), m.colup[j]);
                out << buf;
            }
        }
        if (m.Q.nnz() > 0) {
            out << "QUADOBJ\n";
            for (Idx j = 0; j < m.Q.ncol; ++j) {
                std::string c2 = m.col_names[j];
                for (Off k = m.Q.start[j]; k < m.Q.start[j + 1]; ++k) {
                    Idx i = m.Q.index[k];
                    std::string c1 = m.col_names[i];
                    std::snprintf(buf, sizeof(buf), "    %-10s %-10s %.17g\n", c1.c_str(), c2.c_str(), m.Q.value[k]);
                    out << buf;
                }
            }
        }
        out << "ENDATA\n";
        return true;
    }
};

} // namespace sov

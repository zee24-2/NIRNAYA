NAME          MRPL_Logistics_SupplyChain_MILP
OBJSENSE
  MIN
ROWS
 N  OBJ
 E  ZoneDemand_1
 E  ZoneDemand_2
 E  ZoneDemand_3
 E  ZoneDemand_4
 L  DepotCap_1
 L  DepotCap_2
 L  DepotCap_3
COLUMNS
    MARK0000  'MARKER'                 'INTORG'
    OpenDepot_1 OBJ        1200
    OpenDepot_1 DepotCap_1 -160
    MARK0001  'MARKER'                 'INTEND'
    Ship_d1_z1 OBJ        4
    Ship_d1_z1 ZoneDemand_1 1
    Ship_d1_z1 DepotCap_1 1
    Ship_d1_z2 OBJ        7
    Ship_d1_z2 ZoneDemand_2 1
    Ship_d1_z2 DepotCap_1 1
    Ship_d1_z3 OBJ        8.5
    Ship_d1_z3 ZoneDemand_3 1
    Ship_d1_z3 DepotCap_1 1
    Ship_d1_z4 OBJ        6
    Ship_d1_z4 ZoneDemand_4 1
    Ship_d1_z4 DepotCap_1 1
    MARK0000  'MARKER'                 'INTORG'
    OpenDepot_2 OBJ        900
    OpenDepot_2 DepotCap_2 -140
    MARK0001  'MARKER'                 'INTEND'
    Ship_d2_z1 OBJ        6.5
    Ship_d2_z1 ZoneDemand_1 1
    Ship_d2_z1 DepotCap_2 1
    Ship_d2_z2 OBJ        3.5
    Ship_d2_z2 ZoneDemand_2 1
    Ship_d2_z2 DepotCap_2 1
    Ship_d2_z3 OBJ        5
    Ship_d2_z3 ZoneDemand_3 1
    Ship_d2_z3 DepotCap_2 1
    Ship_d2_z4 OBJ        7.5
    Ship_d2_z4 ZoneDemand_4 1
    Ship_d2_z4 DepotCap_2 1
    MARK0000  'MARKER'                 'INTORG'
    OpenDepot_3 OBJ        1050
    OpenDepot_3 DepotCap_3 -150
    MARK0001  'MARKER'                 'INTEND'
    Ship_d3_z1 OBJ        8
    Ship_d3_z1 ZoneDemand_1 1
    Ship_d3_z1 DepotCap_3 1
    Ship_d3_z2 OBJ        5.5
    Ship_d3_z2 ZoneDemand_2 1
    Ship_d3_z2 DepotCap_3 1
    Ship_d3_z3 OBJ        3
    Ship_d3_z3 ZoneDemand_3 1
    Ship_d3_z3 DepotCap_3 1
    Ship_d3_z4 OBJ        4.5
    Ship_d3_z4 ZoneDemand_4 1
    Ship_d3_z4 DepotCap_3 1
RHS
    RHS1       ZoneDemand_1 70
    RHS1       ZoneDemand_2 85
    RHS1       ZoneDemand_3 65
    RHS1       ZoneDemand_4 60
RANGES
BOUNDS
 BV BND1       OpenDepot_1
 UP BND1       Ship_d1_z1 70
 UP BND1       Ship_d1_z2 85
 UP BND1       Ship_d1_z3 65
 UP BND1       Ship_d1_z4 60
 BV BND1       OpenDepot_2
 UP BND1       Ship_d2_z1 70
 UP BND1       Ship_d2_z2 85
 UP BND1       Ship_d2_z3 65
 UP BND1       Ship_d2_z4 60
 BV BND1       OpenDepot_3
 UP BND1       Ship_d3_z1 70
 UP BND1       Ship_d3_z2 85
 UP BND1       Ship_d3_z3 65
 UP BND1       Ship_d3_z4 60
ENDATA

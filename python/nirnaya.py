"""
NIRNAYA (SIH26119) -- Layer L6 Python Interface
Wraps the native C++17 shared library (bin/nirnaya.dll or bin/libsov.so) via ctypes.
"""
import ctypes
import os
from typing import List, Sequence, Optional


class NirnayaModel:
    CONTINUOUS = 0
    INTEGER = 1
    BINARY = 2

    MINIMIZE = 1
    MAXIMIZE = -1

    def __init__(self, name: str = "nirnaya_model", dll_path: Optional[str] = None):
        if dll_path is None:
            base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            cand = os.path.join(base_dir, "bin", "nirnaya.dll")
            if not os.path.exists(cand):
                cand = os.path.join(base_dir, "bin", "libsov.so")
            dll_path = cand

        self._lib = ctypes.CDLL(dll_path)
        self._configure_signatures()
        self._env = self._lib.sov_env_create()
        self._model = self._lib.sov_model_create(self._env, name.encode("utf-8"))
        self._num_cols = 0

    def _configure_signatures(self):
        self._lib.sov_env_create.restype = ctypes.c_void_p
        self._lib.sov_env_free.argtypes = [ctypes.c_void_p]
        self._lib.sov_model_create.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
        self._lib.sov_model_create.restype = ctypes.c_void_p
        self._lib.sov_model_free.argtypes = [ctypes.c_void_p]

        self._lib.sov_add_col.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_double,
            ctypes.c_double, ctypes.c_int, ctypes.c_char_p
        ]
        self._lib.sov_add_col.restype = ctypes.c_int

        self._lib.sov_add_row.argtypes = [
            ctypes.c_void_p, ctypes.c_double, ctypes.c_double,
            ctypes.c_int, ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_double), ctypes.c_char_p
        ]
        self._lib.sov_add_row.restype = ctypes.c_int

        self._lib.sov_set_sense.argtypes = [ctypes.c_void_p, ctypes.c_int]
        self._lib.sov_set_option_string.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p]
        self._lib.sov_optimize.argtypes = [ctypes.c_void_p]
        self._lib.sov_optimize.restype = ctypes.c_int
        self._lib.sov_get_status.argtypes = [ctypes.c_void_p]
        self._lib.sov_get_status.restype = ctypes.c_int
        self._lib.sov_get_obj.argtypes = [ctypes.c_void_p]
        self._lib.sov_get_obj.restype = ctypes.c_double
        self._lib.sov_get_bound.argtypes = [ctypes.c_void_p]
        self._lib.sov_get_bound.restype = ctypes.c_double
        self._lib.sov_get_x.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double)]

    def add_var(self, cost: float = 0.0, lb: float = 0.0, ub: float = 1e20,
                vtype: int = 0, name: str = "") -> int:
        idx = self._lib.sov_add_col(
            self._model, float(cost), float(lb), float(ub), int(vtype), name.encode("utf-8")
        )
        self._num_cols += 1
        return idx

    def add_constr(self, cols: Sequence[int], vals: Sequence[float],
                   lb: float = -1e20, ub: float = 1e20, name: str = "") -> int:
        nnz = len(cols)
        c_idx = (ctypes.c_int * nnz)(*cols)
        c_val = (ctypes.c_double * nnz)(*vals)
        return self._lib.sov_add_row(
            self._model, float(lb), float(ub), nnz, c_idx, c_val, name.encode("utf-8")
        )

    def set_sense(self, sense: int):
        self._lib.sov_set_sense(self._model, int(sense))

    def set_option(self, key: str, value: str):
        self._lib.sov_set_option_string(self._model, key.encode("utf-8"), str(value).encode("utf-8"))

    def optimize(self) -> int:
        return self._lib.sov_optimize(self._model)

    @property
    def Status(self) -> int:
        return self._lib.sov_get_status(self._model)

    @property
    def ObjVal(self) -> float:
        return self._lib.sov_get_obj(self._model)

    @property
    def ObjBound(self) -> float:
        return self._lib.sov_get_bound(self._model)

    @property
    def X(self) -> List[float]:
        arr = (ctypes.c_double * self._num_cols)()
        self._lib.sov_get_x(self._model, arr)
        return list(arr)

    def __del__(self):
        if getattr(self, "_model", None):
            self._lib.sov_model_free(self._model)
        if getattr(self, "_env", None):
            self._lib.sov_env_free(self._env)

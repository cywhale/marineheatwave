"""No-op compatibility shim (P1-1a).

The MHW API no longer uses the shared Dask scheduler: data is read directly
from zarr with numpy-backed arrays.  These functions are kept so that any
external caller importing them keeps working; they never create a client and
do not import ``distributed``.
"""


def get_dask_client(service_name: str):  # noqa: ARG001
    return None


def close_dask_client(service_name: str):  # noqa: ARG001
    return None

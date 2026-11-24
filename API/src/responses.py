"""Responses that own a temporary file (P1-1b, spec 17).

``FileResponse`` + ``BackgroundTask(os.unlink)`` is not enough: starlette
returns early for a malformed Range (400) or an unsatisfiable Range (416), and
the background task does not run when sending fails, so the file would stay
behind.  ``TempFileResponse`` removes its file in every case once ``__call__``
ends; headers, Range handling and the body come unchanged from FileResponse.
"""
import logging
import os

from starlette.responses import FileResponse

LOG = logging.getLogger(__name__)


def remove_temp_file(path: str, *, during_error: bool = False) -> None:
    """Delete ``path``; a file that is already gone is fine.

    Any other failure is never silently ignored: it is raised, unless another
    exception is already propagating (``during_error``), in which case it is
    logged so that it does not replace the original error.
    """
    try:
        os.unlink(path)
    except FileNotFoundError:
        pass
    except OSError:
        if not during_error:
            raise
        LOG.exception("could not remove temporary file %s", path)


class TempFileResponse(FileResponse):
    """A FileResponse whose file is deleted after the response, however it ends."""

    async def __call__(self, scope, receive, send) -> None:
        try:
            await super().__call__(scope, receive, send)
        except BaseException:
            remove_temp_file(self.path, during_error=True)
            raise
        remove_temp_file(self.path)

"""Low-overhead, flushed diagnostics for verbose network builds."""
import ctypes
import os
import threading
import time
from functools import wraps


def memory_summary():
    """Report native process memory too, without requiring psutil."""
    try:
        if os.name == 'nt':
            from ctypes import wintypes

            class Counters(ctypes.Structure):
                _fields_ = [('cb', wintypes.DWORD), ('faults', wintypes.DWORD)] + [
                    (name, ctypes.c_size_t) for name in (
                        'peak_rss', 'rss', 'peak_paged', 'paged',
                        'peak_nonpaged', 'nonpaged', 'commit', 'peak_commit', 'private')]

            kernel = ctypes.WinDLL('kernel32', use_last_error=True)
            kernel.GetCurrentProcess.restype = wintypes.HANDLE
            query = ctypes.WinDLL('psapi', use_last_error=True).GetProcessMemoryInfo
            query.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
            counters = Counters()
            counters.cb = ctypes.sizeof(counters)
            if not query(kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
                raise OSError(ctypes.get_last_error())
            mib = 1024 ** 2
            return (f'RSS={counters.rss / mib:,.1f} MiB; '
                    f'peak RSS={counters.peak_rss / mib:,.1f} MiB; '
                    f'private={counters.private / mib:,.1f} MiB')
        import resource
        import sys
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        divisor = 1024 ** 2 if sys.platform == 'darwin' else 1024
        return f'peak RSS={peak / divisor:,.1f} MiB'
    except (ImportError, OSError, AttributeError):
        return 'memory unavailable'


def checkpoint(net, stage, details=''):
    if net is None or not net.settings.get('verbose', False):
        return
    net._diagnostic_stage = stage
    net._diagnostic_details = details
    print(f'\n[diagnostic] {stage}' + (f' | {details}' if details else '')
          + f' | {memory_summary()}', flush=True)


def verbose_build(function):
    """Keep reporting the active stage, even inside a slow native operation."""
    @wraps(function)
    def wrapped(net, *args, **kwargs):
        if not net.settings.get('verbose', False):
            return function(net, *args, **kwargs)
        started = time.perf_counter()
        stopped = threading.Event()
        checkpoint(net, 'Network build', f'pid={os.getpid()}; atoms={len(net.balls):,}')

        def heartbeat():
            while not stopped.wait(15):
                print(f'\n[diagnostic] RUNNING {net._diagnostic_stage}'
                      f' | elapsed={time.perf_counter() - started:.1f}s'
                      f' | {net._diagnostic_details} | {memory_summary()}', flush=True)

        worker = threading.Thread(target=heartbeat, daemon=True, name='vorpy-diagnostics')
        worker.start()
        try:
            result = function(net, *args, **kwargs)
        except BaseException as error:
            checkpoint(net, f'FAILED {net._diagnostic_stage}',
                       f'{type(error).__name__}: {error}; {net._diagnostic_details}')
            raise
        else:
            checkpoint(net, 'Network build returned',
                       f'elapsed={time.perf_counter() - started:.1f}s')
            return result
        finally:
            stopped.set()
            worker.join()
    return wrapped

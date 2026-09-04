"""Small first-run installer UI that does not require VorPy's GUI dependencies."""

from __future__ import annotations

import queue
import platform
import subprocess
import sys
import threading
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ICON = Path(__file__).resolve().parent / "assets" / "VorpyIcon_transparent.png"
WINDOWS_ICON = ROOT / "vorpy" / "src" / "GUI" / "Images" / "VorpyIcon.ico"
INSTALL_LOG = ROOT / ".venv" / "vorpy-install.log"
COMMANDS = (
    ([sys.executable, "-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel"],
     "Updating the Python installer..."),
    ([sys.executable, "-m", "pip", "install", "-e", ".[gui]"],
     "Installing VorPy and graphical dependencies..."),
)


def _console_install() -> int:
    for command, message in COMMANDS:
        print(message, flush=True)
        result = subprocess.run(command, cwd=ROOT)
        if result.returncode:
            return result.returncode
    return 0


def _gui_install() -> int:
    import tkinter as tk
    from tkinter import scrolledtext

    events: queue.Queue[tuple[str, object]] = queue.Queue()
    result = {"code": 1}
    if sys.platform == "win32":
        try:
            import ctypes

            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(
                "EricsonLabs.VorPy.Installer"
            )
        except (AttributeError, OSError):
            pass
    root = tk.Tk()
    root.title("VorPy first-time setup")
    root.geometry("680x520")
    root.minsize(520, 400)
    root.configure(bg="#111827")
    if sys.platform == "win32" and WINDOWS_ICON.exists():
        try:
            root.iconbitmap(default=str(WINDOWS_ICON))
        except tk.TclError:
            pass

    try:
        icon = tk.PhotoImage(file=str(ICON))
        root.iconphoto(True, icon)
        scale = max(1, max(icon.width(), icon.height()) // 150)
        logo = icon.subsample(scale, scale)
        logo_label = tk.Label(root, image=logo, bg="#111827")
        logo_label.image = logo
        logo_label.pack(pady=(18, 4))
    except tk.TclError:
        pass

    tk.Label(
        root,
        text="VorPy",
        font=("Helvetica", 22, "bold"),
        fg="#f8fafc",
        bg="#111827",
    ).pack()
    status = tk.StringVar(value="Preparing first-time setup...")
    tk.Label(
        root,
        textvariable=status,
        font=("Helvetica", 11),
        fg="#cbd5e1",
        bg="#111827",
    ).pack(pady=(4, 10))
    log = scrolledtext.ScrolledText(
        root,
        height=12,
        wrap=tk.WORD,
        font=("Menlo", 9),
        bg="#0b1220",
        fg="#dbeafe",
        insertbackground="#dbeafe",
        relief=tk.FLAT,
        padx=10,
        pady=8,
    )
    log.pack(fill=tk.BOTH, expand=True, padx=18, pady=(0, 14))
    close_button = tk.Button(root, text="Close", command=root.destroy, state=tk.DISABLED)
    close_button.pack(pady=(0, 14))
    root.protocol("WM_DELETE_WINDOW", lambda: None)

    def worker() -> None:
        with INSTALL_LOG.open("w", encoding="utf-8") as install_log:
            install_log.write(f"Python: {sys.version}\n")
            install_log.write(f"Platform: {platform.platform()} ({platform.machine()})\n\n")
            for command, message in COMMANDS:
                events.put(("status", message))
                process = subprocess.Popen(
                    command,
                    cwd=ROOT,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                )
                assert process.stdout is not None
                for line in process.stdout:
                    install_log.write(line)
                    install_log.flush()
                    events.put(("line", line))
                code = process.wait()
                if code:
                    events.put(("done", code))
                    return
        events.put(("done", 0))

    def poll() -> None:
        try:
            while True:
                event, value = events.get_nowait()
                if event == "status":
                    status.set(str(value))
                elif event == "line":
                    log.insert(tk.END, str(value))
                    log.see(tk.END)
                elif event == "done":
                    result["code"] = int(value)
                    if value == 0:
                        status.set("Installation complete. Starting VorPy...")
                        root.after(900, root.destroy)
                    else:
                        status.set("Installation failed — review the details below.")
                        close_button.configure(state=tk.NORMAL)
                        root.protocol("WM_DELETE_WINDOW", root.destroy)
                    return
        except queue.Empty:
            pass
        root.after(75, poll)

    threading.Thread(target=worker, daemon=True).start()
    root.after(75, poll)
    root.mainloop()
    return result["code"]


def main() -> int:
    try:
        import tkinter
    except ImportError:
        return _console_install()
    try:
        return _gui_install()
    except tkinter.TclError:
        return _console_install()


if __name__ == "__main__":
    raise SystemExit(main())

import tkinter as tk
from tkinter import ttk
import time


class ProgressWindow(tk.Toplevel):
    def __init__(self, parent):
        super().__init__(parent)
        self.title("Building Network")
        
        # Make window float on top
        self.transient(parent)
        self.grab_set()
        
        # Center the window
        window_width = 400
        window_height = 150
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        x = (screen_width - window_width) // 2
        y = (screen_height - window_height) // 2
        self.geometry(f"{window_width}x{window_height}+{x}+{y}")
        
        # Configure grid
        self.grid_columnconfigure(0, weight=1)
        
        # Create widgets
        self.step_label = ttk.Label(self, text="Initializing...", font=("Arial", 10))
        self.step_label.grid(row=0, column=0, padx=20, pady=(20,5), sticky='w')
        
        self.progress_bar = ttk.Progressbar(self, length=360, mode='determinate')
        self.progress_bar.grid(row=1, column=0, padx=20, pady=5)
        
        stats_frame = ttk.Frame(self)
        stats_frame.grid(row=2, column=0, padx=20, pady=5, sticky='ew')
        stats_frame.grid_columnconfigure(1, weight=1)
        
        ttk.Label(stats_frame, text="Time Elapsed:").grid(row=0, column=0, sticky='w')
        self.time_label = ttk.Label(stats_frame, text="00:00:00")
        self.time_label.grid(row=0, column=1, sticky='e')
        
        ttk.Label(stats_frame, text="Progress:").grid(row=1, column=0, sticky='w')
        self.percent_label = ttk.Label(stats_frame, text="0%")
        self.percent_label.grid(row=1, column=1, sticky='e')
        
        self.start_time = time.time()
        self.update_time()
    
    def update_progress(self, step, progress):
        """Update the progress bar and labels
        
        Args:
            step (str): Current step being performed
            progress (float): Progress percentage (0-100)
        """
        self.step_label.config(text=step)
        self.progress_bar['value'] = progress
        self.percent_label.config(text=f"{progress:.1f}%")
        self.update()
    
    def update_time(self):
        """Update the elapsed time label"""
        elapsed = int(time.time() - self.start_time)
        hours = elapsed // 3600
        minutes = (elapsed % 3600) // 60
        seconds = elapsed % 60
        self.time_label.config(text=f"{hours:02d}:{minutes:02d}:{seconds:02d}")
        self.after(1000, self.update_time)
    
    def finish(self):
        """Close the progress window"""
        self.grab_release()
        self.destroy() 
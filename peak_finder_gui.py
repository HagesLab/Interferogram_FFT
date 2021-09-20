# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 12:33:35 2021

@author: cfai2
"""

import numpy as np
import matplotlib
starting_backend = matplotlib.get_backend()
matplotlib.use("TkAgg")
import matplotlib.backends.backend_tkagg as tkagg
from matplotlib.figure import Figure

import tkinter.filedialog
import tkinter as tk
from tkinter import ttk
import os
from copy import deepcopy as dc
from interferogram_functions import import_MAP
import peak_correction


class tkApp:
    
    def __init__(self, title):
        self.root = tk.Tk()
        self.root.title(title)
        
        s = ttk.Style()
        s.theme_use('classic')
        
        self.data = {"has_data":False, "pos":None, "time":None, "map":None}
        self.peak_params = {"do_BKGsub":True, "BKGlimit":-3,}
        
        self.create_frame()
        
        
    
    def create_frame(self):
        self.main_frame = tk.ttk.Frame(master=self.root).grid(row=0,column=0)
        self.create_loader_frame()
        self.create_plot_frame()
        self.create_sub_controller()
        self.create_subtraction_frame()
        
    def create_loader_frame(self):
        self.loader_frame = tk.ttk.Frame(master=self.main_frame)
        self.loader_frame.grid(row=0,column=0)
        self.do_bkgsub = tk.IntVar()
        self.do_bkgsub.set(1)
        
        tk.ttk.Checkbutton(self.loader_frame, text="Background Subtract", variable=self.do_bkgsub).grid(row=0,column=0,columnspan=2)
        
        tk.ttk.Label(self.loader_frame, text="BKGLimit").grid(row=1,column=0)
        self.BKGlimit_entry = tk.ttk.Entry(self.loader_frame, width=8)
        self.BKGlimit_entry.grid(row=1,column=1)
        self.enter(self.BKGlimit_entry, -3)
        
        tk.ttk.Label(self.loader_frame, text="Search time up to").grid(row=2,column=0)
        self.end_time_tracker = tk.StringVar()
        self.end_time_entry = tk.ttk.Entry(self.loader_frame, textvariable=self.end_time_tracker, width=8)
        
        self.end_time_entry.grid(row=2,column=1)
        self.enter(self.end_time_entry, 20)
        
        tk.ttk.Label(self.loader_frame, text="Search Radius").grid(row=3,column=0)
        self.search_radius_tracker = tk.StringVar()
        self.search_radius_entry = tk.ttk.Entry(self.loader_frame, textvariable=self.search_radius_tracker, width=8)
        
        self.search_radius_entry.grid(row=3, column=1)
        self.enter(self.search_radius_entry, 5)
        
        tk.ttk.Label(self.loader_frame, text="Peak Radius").grid(row=4,column=0)
        self.peak_radius_tracker = tk.StringVar()
        self.peak_radius_entry = tk.ttk.Entry(self.loader_frame, textvariable=self.peak_radius_tracker, width=8)
        
        self.peak_radius_entry.grid(row=4,column=1)
        self.enter(self.peak_radius_entry, 20)
        
        tk.ttk.Label(self.loader_frame, text="Peak Threshold Factor").grid(row=5,column=0)
        self.peak_thr_tracker = tk.StringVar()
        self.peak_thr_entry = tk.ttk.Entry(self.loader_frame, textvariable=self.peak_thr_tracker, width=8)
        
        self.peak_thr_entry.grid(row=5,column=1)
        self.enter(self.peak_thr_entry, 1.5)
        
        tk.ttk.Button(self.loader_frame, text="Load", command=self.load).grid(row=98,column=0,columnspan=2)
        self.update_btn = tk.ttk.Button(self.loader_frame, text="Manual Update", command=self.update_search)
        self.update_btn.grid(row=99,column=0,columnspan=2)

        self.do_auto_update = tk.IntVar()
        self.do_auto_update.trace('w', self.toggle_auto_update)
        self.do_auto_update.set(1)
        tk.ttk.Checkbutton(self.loader_frame, text="Auto update?", variable=self.do_auto_update).grid(row=100,column=0,columnspan=2)
        
    def create_plot_frame(self):
        self.plot_frame = tk.ttk.Frame(master=self.main_frame)
        self.plot_frame.grid(row=0,column=1, rowspan=99)
        
        self.fig = Figure(figsize=(14,8))
        self.search_subplot = self.fig.add_subplot(1,2,1)
        self.search_legend = None
        self.subtract_subplot = self.fig.add_subplot(1,2,2)
        # Prevent coordinate values from appearing in the toolbar; this would sometimes jostle GUI elements around
        #self.search_subplot.format_coord = lambda x, y: ""
        self.canvas = tkagg.FigureCanvasTkAgg(self.fig,master=self.plot_frame)
        self.canvas.get_tk_widget().grid(row=0, column=0)
        self.fig_toolbar_frame = tk.ttk.Frame(master=self.plot_frame)
        self.fig_toolbar_frame.grid(row=1,column=0)
        tkagg.NavigationToolbar2Tk(self.canvas, self.fig_toolbar_frame).grid(row=0,column=0)
        
    def create_subtraction_frame(self):
        self.sub_frame = tk.ttk.Frame(master=self.main_frame)
        self.sub_frame.grid(row=1,column=0)
        
        self.sub_canvas = tk.Canvas(master=self.sub_frame, width=200,height=200)
        self.sub_canvas.grid(row=0,column=0,sticky='nswe')
        
        self.sub_scroll_y = tk.ttk.Scrollbar(self.sub_frame, orient="vertical", 
                                              command=self.sub_canvas.yview)
        self.sub_scroll_y.grid(row=0,column=1, sticky='ns')
        
        self.sub_canvas.configure(yscrollcommand=self.sub_scroll_y.set)
        
        # self.sub_frame.columnconfigure(0,weight=100)
        # self.sub_frame.columnconfigure(1,weight=1, minsize=20) 
        
        self.sub_entry_frame = tk.ttk.Frame(self.sub_canvas)
        
        self.sub_canvas.create_window((0,0), window=self.sub_entry_frame, anchor='nw')
        self.sub_entry_frame.bind('<Configure>', 
                           lambda e:self.sub_canvas.configure(scrollregion=self.sub_canvas.bbox('all')))
        
        
        self.peak_scale_entries = []
        try:
            num_peaks = len(self.peaks) - 1
        except Exception:
            num_peaks = 0
        for n in range(num_peaks):
            tk.ttk.Label(self.sub_entry_frame, text="Peak #{}".format(n+1)).grid(row=n,column=0)
            self.peak_scale_entries.append(tk.ttk.Entry(self.sub_entry_frame, width=8))
            self.peak_scale_entries[-1].grid(row=n,column=1)
            
        self.toggle_auto_subtract()
                
    def create_sub_controller(self):
        self.sub_control_frame = tk.ttk.Frame(master=self.main_frame)
        self.sub_control_frame.grid(row=2,column=0)
        self.subtract_btn = tk.ttk.Button(self.sub_control_frame, text="Subtract", command=self.subtract_peaks)
        self.subtract_btn.grid(row=0,column=0,columnspan=2)
        
        self.do_auto_subtract = tk.IntVar()
        self.do_auto_subtract.set(1)
        self.do_auto_subtract.trace('w', self.toggle_auto_subtract)
        tk.ttk.Checkbutton(self.sub_control_frame, text="Auto subtract?", variable=self.do_auto_subtract).grid(row=1,column=0,columnspan=2)
        
        tk.ttk.Button(self.sub_control_frame, text="Save", command=self.save).grid(row=2,column=0,columnspan=2)
        
        
    def shift_time_to_max(self):
        t_max = self.data["time"][np.array(np.where(np.mean(self.data["map"],axis=0)==np.max(np.mean(self.data["map"],axis=0)))[0],dtype="int")]
        self.data["time"] -= t_max
        return
        
    def remove_background(self):
        # TODO: Popup result of this
        self.peak_params["do_BKGsub"] = self.do_bkgsub.get()
        self.peak_params["BKGlimit"] = float(self.BKGlimit_entry.get())
        
        if self.peak_params["do_BKGsub"]:
            BKGrange = np.array([self.data["time"][0], self.peak_params["BKGlimit"]],dtype='float')  #ns
            
            index = [(np.abs(self.data["time"]-np.min(BKGrange))).argmin(),(np.abs(self.data["time"]-np.max(BKGrange))).argmin()]
            BKGval = np.mean(self.data["map"][:,np.min(index):np.max(index)],axis=1)
            self.data["map"] -= np.reshape(BKGval, (len(BKGval), 1))
            
        
    def load(self):
        print("Loaded plot")
        self.data["has_data"] = False
        try:        
            fname = tk.filedialog.askdirectory(title="Select MAP data directory")
            self.data["pos"], self.data["time"], self.data["map"] = import_MAP(fname)
            self.data["has_data"] = True
        except Exception:
            print("Error: could not read {}".format(fname))
            return
        
        self.shift_time_to_max()
        self.remove_background()
        self.update_search()
        self.subtract_peaks()
        return
    
    def toggle_auto_update(self, *args):
        if self.do_auto_update.get():
            self.update_btn.configure(state='disabled')
            self.end_time_entry.bind("<Return>", self.update_search)
            self.search_radius_entry.bind("<Return>", self.update_search)
            self.peak_radius_entry.bind("<Return>", self.update_search)
            self.peak_thr_entry.bind("<Return>", self.update_search)
        else:
            self.update_btn.configure(state='enabled')
            self.end_time_entry.unbind("<Return>")
            self.search_radius_entry.unbind("<Return>")
            self.peak_radius_entry.unbind("<Return>")
            self.peak_thr_entry.unbind("<Return>")
        
    def update_search(self, *args):
        # TODO: Have this happen whenever a peak_params is changed
        if self.data["has_data"]:
            self.update_search_plot()
            self.sub_frame.grid_forget()
            self.create_subtraction_frame()
        else:
            print("No data yet")
        
        return
    
    def update_search_plot(self):
                 
        TRPLmin_OM = 1e-6
        try:
            self.peak_params["end_time"] = float(self.end_time_tracker.get())
        except Exception:
            print("Invalid search time")
            return
        
        try:
            self.peak_params["search_radius"] = int(self.search_radius_tracker.get())
        except Exception:
            print("Invalid search radius")
            return
        
        try:
            self.peak_params["peak_radius"] = int(self.peak_radius_tracker.get())
        except Exception:
            print("Invalid peak radius")
            return
        
        try:
            self.peak_params["peak_thr"] = float(self.peak_thr_tracker.get())
        except Exception:
            print("Invalid peak threshold")
            return
        
        # Optional: set to FALSE to hide this plot when no longer needed
        plot_BKGsub = True
        
    
        integralTRPL = np.sum(self.data["map"],axis=0)
        # Ensure that the arguments of find_peaks are finding the peaks correctly first
        self.peaks = peak_correction.find_peaks(integralTRPL, self.data["time"], end_time=self.peak_params["end_time"],
                                           search_radius=self.peak_params["search_radius"], peak_thr=self.peak_params["peak_thr"],
                                           peak_radius=self.peak_params["peak_radius"])
            
        #Plot Full TRPL
        self.search_subplot.cla()
        if self.search_legend is not None: self.search_legend.remove()
        self.search_subplot.set_title("Peaks Identified")
        self.search_subplot.plot(self.data["time"],integralTRPL)
        self.search_subplot.set_ylim(np.max(integralTRPL)*TRPLmin_OM,2*np.max(integralTRPL))
        self.search_subplot.set_xlabel('Time / ns')
        self.search_subplot.set_ylabel('Counts / a.u.')
        self.search_subplot.set_yscale('log')
        
        print("Found peaks: ", self.peaks)
        for i, peak in enumerate(self.peaks):
            if i > 0:
                label= "Peak #{}".format(i)
            else:
                label = "Main Peak"
            self.search_subplot.plot(self.data["time"][peak[0]:peak[1]], integralTRPL[peak[0]:peak[1]], label=label)
            
        self.fig.tight_layout()
        self.search_legend = self.fig.legend(loc='upper left')
        self.search_legend.set_draggable(True)
        self.fig.canvas.draw()
        
    def toggle_auto_subtract(self, *args):
        if self.do_auto_subtract.get():
            self.subtract_btn.configure(state='disabled')
            for entrybox in self.peak_scale_entries:
                entrybox.bind("<Return>", self.subtract_peaks)
        else:
            self.subtract_btn.configure(state='enabled')
            for entrybox in self.peak_scale_entries:
                entrybox.unbind("<Return>")
    
    def subtract_peaks(self, *args):
        if self.data["has_data"]:
            integralTRPL = np.sum(self.data["map"],axis=0)
            TRPLmin_OM = 1e-6
            
            self.reduce_factors = []
            for n in self.peak_scale_entries:
                try:
                    self.reduce_factors.append(float(n.get()))
                except Exception:
                    print("Invalid peak scale {}; skipping".format(n.get()))
                    self.reduce_factors.append(0)
                    
            peak_correction.subtract_peaks(integralTRPL, self.peaks, reduce=self.reduce_factors)
            
            self.subtract_subplot.cla()
            self.subtract_subplot.set_title("Peaks Subtracted")
            self.subtract_subplot.plot(self.data["time"],integralTRPL)
            self.subtract_subplot.set_ylim(np.max(integralTRPL)*TRPLmin_OM,2*np.max(integralTRPL))
            self.subtract_subplot.set_xlabel('Time / ns')
            self.subtract_subplot.set_ylabel('Counts / a.u.')
            self.subtract_subplot.set_yscale('log')
            
            self.fig.tight_layout()
            self.fig.canvas.draw()
            
    def save(self):
        if self.data["has_data"]:
            new_map = np.array(self.data['map'])
            peaks_temp = dc(self.peaks)
            for wavelength in new_map:
                # This overwrites self.peaks
                peak_correction.find_peaks(wavelength, self.data['time'], peaks_known=self.peaks)
                peak_correction.subtract_peaks(wavelength, self.peaks, reduce=self.reduce_factors)
            
            path = tk.filedialog.asksaveasfilename(title="Save corrected map", 
                                                   filetypes=[("Text files","*.txt")])
            
            
            if path:
                if path.endswith(".txt"): 
                    path = path[:-4]
                np.savetxt("{}.txt".format(path), new_map, delimiter="\t")
        
            # Regenerate original self.peaks
            self.peaks = dc(peaks_temp)
            
    def run(self):
        # width, height = self.root.winfo_screenwidth() * 0.8, self.root.winfo_screenheight() * 0.8

        # self.root.geometry('%dx%d+0+0' % (width,height))
        self.root.attributes("-topmost", True)
        self.root.after_idle(self.root.attributes,'-topmost',False)
        self.root.mainloop()
        print("Closed")
        matplotlib.use(starting_backend)
        return

    def quit(self):
        self.root.destroy()
        print("Closed")
        matplotlib.use(starting_backend)
        return
    
    def enter(self, entryBox, text):
        """ Fill user entry boxes with text. """
        entryBox.delete(0,tk.END)
        entryBox.insert(0,text)
        return

if __name__ == "__main__":
    tkapp = tkApp("Peak Finder")
    tkapp.run()
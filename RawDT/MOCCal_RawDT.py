# MOCCal Script for Biomolecular Class Assignment and CCS Calibration

from tkinter import *
from tkinter import ttk
from tkinter.filedialog import askopenfile 
from tkinter.filedialog import askdirectory
from tkinter import simpledialog
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import math

#EDC GUI

ws = Tk()
ws.withdraw()

edc = simpledialog.askfloat("EDC", "What is your EDC delay coefficient?")

ws.destroy()

ws.mainloop()

#Class Selection + Calibration Solution Selection + Data Import GUI

ws = Tk()
ws.title('MOCCal Data Selection')
ws.geometry('600x500') 
 
check_1 = IntVar()
check_2 = IntVar()
check_3 = IntVar()
 
lipids_present = 0
metabolites_present = 0
peptides_present = 0
 
def lipidbtn_Clicked():
    global lipids_present
    if check_1.get() :
        lipids_present = 1
    else :
        lipids_present = 0
 
def metabbtn_Clicked():
    global metabolites_present
    if check_2.get() :
        metabolites_present = 1
    else :
        metabolites_present = 0
 
def pepbtn_Clicked():
    global peptides_present
    if check_3.get() :
        peptides_present = 1
    else :
        peptides_present = 0

def nextPage():
    ws.destroy()
        
Label(ws, text='Please Select Classes Present', foreground='green').grid(row=0, column=1, padx=10)
 
lipid_btn = Checkbutton(ws, text = 'Lipids', variable = check_1, onvalue = 1, offvalue = 0, height=2, width = 20, command = lipidbtn_Clicked)
lipid_btn.grid(row=1, column=0, padx=10, pady=10)
 
metab_btn = Checkbutton(ws, text = 'Small Molecules', variable = check_2, onvalue = 1, offvalue = 0, height=2, width = 20, command = metabbtn_Clicked)
metab_btn.grid(row=1, column=1, padx=10, pady=10)
 
pep_btn = Checkbutton(ws, text = 'Peptides', variable = check_3, onvalue = 1, offvalue = 0, height=2, width = 20, command = pepbtn_Clicked)
pep_btn.grid(row=1, column=2, padx=10, pady=10)
 
solution_check_1 = IntVar()
solution_check_2 = IntVar()
 
combined_present = 0
individual_present = 0
 
def combobtn_Clicked():
    global combined_present
    if solution_check_1.get() :
        combined_present = 1
    else :
        combined_present = 0
 
def indivbtn_Clicked():
    global individual_present
    if solution_check_2.get() :
        individual_present = 1
    else :
        individual_present = 0
        
Label(ws, text='Please Select Calibration Solution Option', foreground='green').place(relx=0.3, rely=0.2, anchor=W)
 
combined_btn = Checkbutton(ws, text = 'Calibrants in Combined Solution', variable = solution_check_1, onvalue = 1, offvalue = 0, height=2, width = 30, command = combobtn_Clicked)
combined_btn.place(relx = 0.09, rely = 0.29, anchor= W)
 
individual_btn = Checkbutton(ws, text = 'Calibrants in Individual Solutions', variable = solution_check_2, onvalue = 1, offvalue = 0, height=2, width = 30, command = indivbtn_Clicked)
individual_btn.place(relx = 0.5, rely = 0.29, anchor= W)

Label(ws, text='Please Select Calibration Excel Sheet', foreground='green').place(relx = 0.32, rely = 0.4, anchor= W)

def on_calibrant_excel_file_click():
    file_path = askopenfile(mode='r', filetypes=[('Excel Files', '*.xlsx')])
    if file_path is not None:
        global calibrant_file_name
        calibrant_file_name = str(file_path.name)
        calibrant_label["text"] = calibrant_file_name

calibrant_label = Label(ws, text='Calibration Excel Sheet')
calibrant_label.place(relx = 0.23, rely = 0.49, anchor= W)

calibrant_btn = Button(ws, text ='Choose File', 
    command = lambda:on_calibrant_excel_file_click()) 
calibrant_btn.place(relx = 0.1, rely = 0.49, anchor= W)

def on_experimental_excel_file_click():
    file_path = askopenfile(mode='r', filetypes = [('Excel Files', ('*.xlsx', '*.xls'))])
    if file_path is not None:
        global exp_data_name
        exp_data_name = str(file_path.name)
        exp_data_label["text"] = exp_data_name

Label(ws, text='Please Select Experimental Data', foreground='green').place(relx=0.345, rely=0.6, anchor=W)
 
exp_data_label = Label(ws, text='Experimental Data Excel Sheet')
exp_data_label.place(relx = 0.23, rely = 0.69, anchor= W)

exp_data_btn = Button(ws, text ='Choose File', 
    command = lambda:on_experimental_excel_file_click()) 
exp_data_btn.place(relx = 0.1, rely = 0.69, anchor= W)

next_btn = Button(ws, text ='Next', command = nextPage).place(relx = 0.5, rely = 0.85, anchor= S)

ws.mainloop()

if individual_present == 1:

    ws = Tk()
    ws.title('Import Data')
    ws.geometry('700x200') 
    
    def on_lipid_calibration_data_file_click():
       file_path = askopenfile(mode='r', filetypes=[('mzML Files', '*.mzML')])
       if file_path is not None:      
            global lipid_data_file
            lipid_data_file = str(file_path.name)
            lipid_calibration_data_label["text"] = lipid_data_file
    def on_metabolite_calibration_data_file_click():
       file_path = askopenfile(mode='r', filetypes=[('mzML Files', '*.mzML')])
       if file_path is not None:       
            global metabolite_data_file
            metabolite_data_file = str(file_path.name)
            metabolite_calibration_data_label["text"] = metabolite_data_file
    def on_peptide_calibration_data_file_click():
       file_path = askopenfile(mode='r', filetypes=[('mzML Files', '*.mzML')])
       if file_path is not None:
            global peptide_data_file
            peptide_data_file = str(file_path.name)
            peptide_calibration_data_label["text"] = peptide_data_file
            
    if lipids_present == 1:
    
        lipid_calibration_data_label = Label(ws, text='Upload Lipid Calibration Data as .mzML')
        lipid_calibration_data_label.place(relx = 0.2, rely = 0.2, anchor= W)
    
        lipid_calibration_data_btn = Button(ws, text ='Choose File', 
            command = lambda:on_lipid_calibration_data_file_click()) 
        lipid_calibration_data_btn.place(relx = 0.05, rely = 0.2, anchor= W)
              
    if metabolites_present == 1:
    
        metabolite_calibration_data_label = Label(ws, text='Upload Small Molecule Calibration Data as .mzML')
        if lipids_present == 1:
            metabolite_calibration_data_label.place(relx = 0.2, rely = 0.4, anchor= W)
        else: metabolite_calibration_data_label.place(relx = 0.2, rely = 0.2, anchor= W)
    
        metabolite_calibration_data_btn = Button(ws, text ='Choose File', 
            command = lambda:on_metabolite_calibration_data_file_click())
        if lipids_present == 1:
            metabolite_calibration_data_btn.place(relx = 0.05, rely = 0.4, anchor= W)
        else: metabolite_calibration_data_btn.place(relx = 0.05, rely = 0.2, anchor= W)
        
    if peptides_present == 1:
    
        peptide_calibration_data_label = Label(ws, text='Upload Peptide Calibration Data as .mzML')
        if all([lipids_present, metabolites_present]):
            peptide_calibration_data_label.place(relx = 0.2, rely = 0.6, anchor= W)
        elif any([lipids_present, metabolites_present]) and (all([lipids_present, metabolites_present]) != True):
            peptide_calibration_data_label.place(relx = 0.2, rely = 0.4, anchor= W)
        else: peptide_calibration_data_label.place(relx = 0.2, rely = 0.2, anchor= W)
    
        peptide_calibration_data_btn = Button(ws, text ='Choose File', 
            command = lambda:on_peptide_calibration_data_file_click()) 
        if all([lipids_present, metabolites_present]):
            peptide_calibration_data_btn.place(relx = 0.05, rely = 0.6, anchor= W)
        elif any([lipids_present, metabolites_present]) and (all([lipids_present, metabolites_present]) != True):
            peptide_calibration_data_btn.place(relx = 0.05, rely = 0.4, anchor= W)
        else: peptide_calibration_data_btn.place(relx = 0.05, rely = 0.2, anchor= W)
     
    next_btn = Button(ws, text ='Process', command = nextPage).place(relx = 0.5, rely = 1.0, anchor= S)

    ws.mainloop()

elif combined_present == 1:
    ws = Tk()
    ws.title('Import Data')
    ws.geometry('700x100') 
    
    def on_calibration_data_file_click():
       file_path = askopenfile(mode='r', filetypes=[('mzML Files', '*.mzML')])
       if file_path is not None:      
            global calibrant_data_file
            calibrant_data_file = str(file_path.name)
            calibration_data_label["text"] = calibrant_data_file
    
    calibration_data_label = Label(ws, text='Upload Calibration Data as .mzML')
    calibration_data_label.place(relx = 0.2, rely = 0.2, anchor= W)

    calibration_data_btn = Button(ws, text ='Choose File', 
        command = lambda:on_calibration_data_file_click()) 
    calibration_data_btn.place(relx = 0.05, rely = 0.2, anchor= W)
    next_btn = Button(ws, text ='Process', command = nextPage).place(relx = 0.5, rely = 1.0, anchor= S)
    
    ws.mainloop()
    
all_calibrants = pd.read_excel(calibrant_file_name)
lipid_calibrants = all_calibrants[all_calibrants["biomolecule"].isin(["lipid"])]
metabolite_calibrants = all_calibrants[all_calibrants["biomolecule"].isin(["small molecule"])]
pep_calibrants_all = all_calibrants[all_calibrants["biomolecule"].isin(["peptide"])]  

if combined_present == 1:
    lipid_data_file = calibrant_data_file
    metabolite_data_file = calibrant_data_file
    peptide_data_file = calibrant_data_file
    
import deimos
import numpy as np
from scipy.optimize import curve_fit
import math

nitrogen_mass = 28.0134

cal_parameter_names = ['A', 't0', 'B']

lipid_A = 0
if lipids_present == 1:

    calibrants = lipid_calibrants
    
    accessions = deimos.get_accessions(lipid_data_file)
    
    data = deimos.load(lipid_data_file,
                       accession={'retention_time': 'MS:1000016',
                                  'drift_time': 'MS:1002476'})
    
    ms1 = data['ms1']
    ms1= ms1.drop(['function', 'process', 'scan'], axis=1)
    
    #nominal threshold
    ms1 = deimos.threshold(ms1, threshold=128)
    
    #sum over retention time
    ms1_sum = deimos.collapse(ms1, keep=['mz', 'drift_time'])
    ms1_sum= ms1_sum.sort_values(by='mz', ascending=True)
    
    calibrant_name = np.array(calibrants['calibrant'])
    mass = np.array(calibrants['m/z'])
    lit_ccs = np.array(calibrants['ccs'])
    
    extracted_drift_time = np.empty(shape= [0,1])
    
    def Gauss(x, amplitude, mean, stdev):
        y = amplitude * np.exp(-0.5*((x-mean)/stdev)**2)
        return y
    
    for i in range(0,len(mass)):
        ms_slice = deimos.slice(ms1_sum, by='mz', low=[mass[i]-0.05], high= [mass[i]+0.05])
        slice_peaks = deimos.peakpick.local_maxima(ms_slice, dims=['mz', 'drift_time'], bins=[37, 1])
        slice_peaks = slice_peaks.sort_values(by='drift_time', ascending=True)
        dt = np.array(slice_peaks['drift_time'])
        intensity = np.array(slice_peaks['intensity']) 
        parameters, covariance = curve_fit(Gauss, dt, intensity, maxfev=1000000, p0 = (max(intensity), dt[np.argmax(intensity)], 0.5))
        Gauss_intensity = round(parameters[0],0)
        Gauss_dt = round(parameters[1], 2)
        Gauss_stdev = round(parameters[2],2)
        extracted_drift_time = np.append(extracted_drift_time, Gauss_dt)
    
    # Calculating CCS Calibration Parameters
    
    corrected_dt = []
    corrected_ccs = []
    reduced_mass = []
    for i in range(0, len(mass)):
        reduced_mass_val = (mass[i] * nitrogen_mass)/(mass[i] + nitrogen_mass)
        reduced_mass = np.append(reduced_mass, reduced_mass_val)
        corrected_dt_val = extracted_drift_time[i] - (((math.sqrt(mass[i]))*edc)/1000)
        corrected_ccs_val = lit_ccs[i] / math.sqrt(1/reduced_mass_val)
        corrected_dt = np.append(corrected_dt, corrected_dt_val)
        corrected_ccs = np.append(corrected_ccs, corrected_ccs_val)
    
    def cal_curve(corrected_dt, A, t0, B):
        y = A * (corrected_dt + t0)**B
        return y
    
    param, cov = curve_fit(cal_curve, corrected_dt, corrected_ccs, maxfev=1000000, p0=(500., 0.001, 0.5))
    
    lipid_A = param[0]
    lipid_t0 = param[1]
    lipid_B = param[2]
    
    # CALCULATES CCS VALUES AND RESIDUALS
    
    cal_ccs = []
    for i in range (0, len(corrected_dt)):
        cal_ccs_val = lipid_A * math.sqrt(1/reduced_mass[i]) * (corrected_dt[i] + lipid_t0)**lipid_B
        cal_ccs = np.append(cal_ccs, round(cal_ccs_val,2))
    
    residual = []
    for i in range(len(lit_ccs)):
        new_residual = round(((cal_ccs[i]-lit_ccs[i])/lit_ccs[i])*100,3)
        residual = np.append(residual, new_residual)
        
    # FINAL DATA OUTPUT
    
    lipid_final_data = pd.DataFrame({'calibrant': calibrant_name, 'mz': mass, 'Drift Time': extracted_drift_time, 'Lit CCS': lit_ccs, 'Cal CCS': cal_ccs, 'Residuals': residual})
    lipid_parameters = pd.DataFrame({'parameter': cal_parameter_names, 'values': [lipid_A, lipid_t0, lipid_B]})
    

metab_A = 0
if metabolites_present == 1: 
    
    calibrants = metabolite_calibrants
    
    accessions = deimos.get_accessions(metabolite_data_file)
    
    data = deimos.load(metabolite_data_file,
                       accession={'retention_time': 'MS:1000016',
                                  'drift_time': 'MS:1002476'})
    
    ms1 = data['ms1']
    ms1= ms1.drop(['function', 'process', 'scan'], axis=1)
    
    #nominal threshold
    ms1 = deimos.threshold(ms1, threshold=128)
    
    #sum over retention time
    ms1_sum = deimos.collapse(ms1, keep=['mz', 'drift_time'])
    ms1_sum= ms1_sum.sort_values(by='mz', ascending=True)
    
    calibrant_name = np.array(calibrants['calibrant'])
    mass = np.array(calibrants['m/z'])
    lit_ccs = np.array(calibrants['ccs'])
    
    extracted_drift_time = np.empty(shape= [0,1])
    
    def Gauss(x, amplitude, mean, stdev):
        y = amplitude * np.exp(-0.5*((x-mean)/stdev)**2)
        return y
    
    for i in range(0,len(mass)):
        ms_slice = deimos.slice(ms1_sum, by='mz', low=[mass[i]-0.05], high= [mass[i]+0.05])
        slice_peaks = deimos.peakpick.local_maxima(ms_slice, dims=['mz', 'drift_time'], bins=[37, 1])
        slice_peaks = slice_peaks.sort_values(by='drift_time', ascending=True)
        dt = np.array(slice_peaks['drift_time'])
        intensity = np.array(slice_peaks['intensity']) 
        parameters, covariance = curve_fit(Gauss, dt, intensity, maxfev=1000000, p0 = (max(intensity), dt[np.argmax(intensity)], 0.5))
        Gauss_intensity = round(parameters[0],0)
        Gauss_dt = round(parameters[1], 2)
        Gauss_stdev = round(parameters[2],2)
        extracted_drift_time = np.append(extracted_drift_time, Gauss_dt)
    
    # Calculating CCS Calibration Parameters
    
    corrected_dt = []
    corrected_ccs = []
    reduced_mass = []
    for i in range(0, len(mass)):
        reduced_mass_val = (mass[i] * nitrogen_mass)/(mass[i] + nitrogen_mass)
        reduced_mass = np.append(reduced_mass, reduced_mass_val)
        corrected_dt_val = extracted_drift_time[i] - (((math.sqrt(mass[i]))*edc)/1000)
        corrected_ccs_val = lit_ccs[i] / math.sqrt(1/reduced_mass_val)
        corrected_dt = np.append(corrected_dt, corrected_dt_val)
        corrected_ccs = np.append(corrected_ccs, corrected_ccs_val)
    
    def cal_curve(corrected_dt, A, t0, B):
        y = A * (corrected_dt + t0)**B
        return y
    
    param, cov = curve_fit(cal_curve, corrected_dt, corrected_ccs, maxfev=1000000, p0=(500., 0.001, 0.5))
    
    metab_A = param[0]
    metab_t0 = param[1]
    metab_B = param[2]
    
    # CALCULATES CCS VALUES AND RESIDUALS
    
    cal_ccs = []
    for i in range (0, len(corrected_dt)):
        cal_ccs_val = metab_A * math.sqrt(1/reduced_mass[i]) * (corrected_dt[i] + metab_t0)**metab_B
        cal_ccs = np.append(cal_ccs, round(cal_ccs_val,2))
    
    residual = []
    for i in range(len(lit_ccs)):
        new_residual = round(((cal_ccs[i]-lit_ccs[i])/lit_ccs[i])*100,3)
        residual = np.append(residual, new_residual)
        
    # FINAL DATA OUTPUT
    
    metabolite_final_data = pd.DataFrame({'calibrant': calibrant_name, 'mz': mass, 'Drift Time': extracted_drift_time, 'Lit CCS': lit_ccs, 'Cal CCS': cal_ccs, 'Residuals': residual})
    metabolite_parameters = pd.DataFrame({'parameter': cal_parameter_names, 'values': [metab_A, metab_t0, metab_B]})

single_pep_A = 0
double_pep_A = 0
triple_pep_A = 0
if peptides_present == 1: 
    
    single_data = pep_calibrants_all.loc[pep_calibrants_all['charge'] == 1]
    double_data = pep_calibrants_all.loc[pep_calibrants_all['charge'] == 2]
    triple_data = pep_calibrants_all.loc[pep_calibrants_all['charge'] == 3]
    
    single_pep_A = 0
    double_pep_A = 0
    triple_pep_A = 0
    
    if any(pep_calibrants_all['charge']==1):

        calibrants = single_data
        accessions = deimos.get_accessions(peptide_data_file)
        data = deimos.load(peptide_data_file,
                           accession={'retention_time': 'MS:1000016',
                                      'drift_time': 'MS:1002476'})
        ms1 = data['ms1']
        ms1= ms1.drop(['function', 'process', 'scan'], axis=1)
        
        #nominal threshold
        ms1 = deimos.threshold(ms1, threshold=128)
        
        #sum over retention time
        ms1_sum = deimos.collapse(ms1, keep=['mz', 'drift_time'])
        ms1_sum= ms1_sum.sort_values(by='mz', ascending=True)
        
        calibrant_name = np.array(calibrants['calibrant'])
        mass = np.array(calibrants['m/z'])
        lit_ccs = np.array(calibrants['ccs'])
        
        extracted_drift_time = np.empty(shape= [0,1])
        
        def Gauss(x, amplitude, mean, stdev):
            y = amplitude * np.exp(-0.5*((x-mean)/stdev)**2)
            return y
        
        for i in range(0,len(mass)):
            ms_slice = deimos.slice(ms1_sum, by='mz', low=[mass[i]-0.05], high= [mass[i]+0.05])
            slice_peaks = deimos.peakpick.local_maxima(ms_slice, dims=['mz', 'drift_time'], bins=[37, 1])
            slice_peaks = slice_peaks.sort_values(by='drift_time', ascending=True)
            dt = np.array(slice_peaks['drift_time'])
            intensity = np.array(slice_peaks['intensity']) 
            parameters, covariance = curve_fit(Gauss, dt, intensity, maxfev=1000000, p0 = (max(intensity), dt[np.argmax(intensity)], 0.5))
            Gauss_intensity = round(parameters[0],0)
            Gauss_dt = round(parameters[1], 2)
            Gauss_stdev = round(parameters[2],2)
            extracted_drift_time = np.append(extracted_drift_time, Gauss_dt)
        
        # Calculating CCS Calibration Parameters
        
        corrected_dt = []
        corrected_ccs = []
        reduced_mass = []
        for i in range(0, len(mass)):
            reduced_mass_val = (mass[i] * nitrogen_mass)/(mass[i] + nitrogen_mass)
            reduced_mass = np.append(reduced_mass, reduced_mass_val)
            corrected_dt_val = extracted_drift_time[i] - (((math.sqrt(mass[i]))*edc)/1000)
            corrected_ccs_val = lit_ccs[i] / math.sqrt(1/reduced_mass_val)
            corrected_dt = np.append(corrected_dt, corrected_dt_val)
            corrected_ccs = np.append(corrected_ccs, corrected_ccs_val)
        
        def cal_curve(corrected_dt, A, t0, B):
            y = A * (corrected_dt + t0)**B
            return y
        
        param, cov = curve_fit(cal_curve, corrected_dt, corrected_ccs, maxfev=1000000, p0=(500., 0.001, 0.5))
        
        single_pep_A = param[0]
        single_pep_t0 = param[1]
        single_pep_B = param[2] 
        
        # CALCULATES CCS VALUES AND RESIDUALS
        
        cal_ccs = []
        for i in range (0, len(corrected_dt)):
            cal_ccs_val = (single_pep_A * math.sqrt(1/reduced_mass[i]) * (corrected_dt[i] + single_pep_t0)**single_pep_B)*1
            cal_ccs = np.append(cal_ccs, round(cal_ccs_val,2))
        
        residual = []
        for i in range(len(lit_ccs)):
            new_residual = round(((cal_ccs[i]-lit_ccs[i])/lit_ccs[i])*100,3)
            residual = np.append(residual, new_residual)
            
        # FINAL DATA OUTPUT
        
        single_peptide_final_data = pd.DataFrame({'calibrant': calibrant_name, 'mz': mass, 'Drift Time': extracted_drift_time, 'Lit CCS': lit_ccs, 'Cal CCS': cal_ccs, 'Residuals': residual})
        single_peptide_parameters = pd.DataFrame({'parameter': cal_parameter_names, 'values': [single_pep_A, single_pep_t0, single_pep_B]})
        
    if any(pep_calibrants_all['charge']==2):

        calibrants = double_data
        accessions = deimos.get_accessions(peptide_data_file)
        data = deimos.load(peptide_data_file,
                           accession={'retention_time': 'MS:1000016',
                                      'drift_time': 'MS:1002476'})
        ms1 = data['ms1']
        ms1= ms1.drop(['function', 'process', 'scan'], axis=1)
        
        #nominal threshold
        ms1 = deimos.threshold(ms1, threshold=128)
        
        #sum over retention time
        ms1_sum = deimos.collapse(ms1, keep=['mz', 'drift_time'])
        ms1_sum= ms1_sum.sort_values(by='mz', ascending=True)
        
        calibrant_name = np.array(calibrants['calibrant'])
        mass = np.array(calibrants['m/z'])
        lit_ccs = np.array(calibrants['ccs'])
        
        extracted_drift_time = np.empty(shape= [0,1])
        
        def Gauss(x, amplitude, mean, stdev):
            y = amplitude * np.exp(-0.5*((x-mean)/stdev)**2)
            return y
        
        for i in range(0,len(mass)):
            ms_slice = deimos.slice(ms1_sum, by='mz', low=[mass[i]-0.05], high= [mass[i]+0.05])
            slice_peaks = deimos.peakpick.local_maxima(ms_slice, dims=['mz', 'drift_time'], bins=[37, 1])
            slice_peaks = slice_peaks.sort_values(by='drift_time', ascending=True)
            dt = np.array(slice_peaks['drift_time'])
            intensity = np.array(slice_peaks['intensity']) 
            parameters, covariance = curve_fit(Gauss, dt, intensity, maxfev=1000000, p0 = (max(intensity), dt[np.argmax(intensity)], 0.5))
            Gauss_intensity = round(parameters[0],0)
            Gauss_dt = round(parameters[1], 2)
            Gauss_stdev = round(parameters[2],2)
            extracted_drift_time = np.append(extracted_drift_time, Gauss_dt)
        
        # Calculating CCS Calibration Parameters
        
        corrected_dt = []
        corrected_ccs = []
        reduced_mass = []
        for i in range(0, len(mass)):
            reduced_mass_val = (mass[i] * nitrogen_mass)/(mass[i] + nitrogen_mass)
            reduced_mass = np.append(reduced_mass, reduced_mass_val)
            corrected_dt_val = extracted_drift_time[i] - (((math.sqrt(mass[i]))*edc)/1000)
            corrected_ccs_val = lit_ccs[i] * (math.sqrt(reduced_mass_val))/2
            corrected_dt = np.append(corrected_dt, corrected_dt_val)
            corrected_ccs = np.append(corrected_ccs, corrected_ccs_val)
        
        def cal_curve(corrected_dt, A, t0, B):
            y = A * (corrected_dt + t0)**B
            return y
        
        param, cov = curve_fit(cal_curve, corrected_dt, corrected_ccs, maxfev=1000000, p0=(500., 0.001, 0.5))
        
        double_pep_A = param[0]
        double_pep_t0 = param[1]
        double_pep_B = param[2]
        
        # CALCULATES CCS VALUES AND RESIDUALS
        
        cal_ccs = []
        for i in range (0, len(corrected_dt)):
            cal_ccs_val = (double_pep_A * math.sqrt(1/reduced_mass[i]) * (corrected_dt[i] + double_pep_t0)**double_pep_B)*2
            cal_ccs = np.append(cal_ccs, round(cal_ccs_val,2))
        
        residual = []
        for i in range(len(lit_ccs)):
            new_residual = round(((cal_ccs[i]-lit_ccs[i])/lit_ccs[i])*100,3)
            residual = np.append(residual, new_residual)
            
        # FINAL DATA OUTPUT
        
        double_peptide_final_data = pd.DataFrame({'calibrant': calibrant_name, 'mz': mass, 'Drift Time': extracted_drift_time, 'Lit CCS': lit_ccs, 'Cal CCS': cal_ccs, 'Residuals': residual})
        double_peptide_parameters = pd.DataFrame({'parameter': cal_parameter_names, 'values': [double_pep_A, double_pep_t0, double_pep_B]})
        #print(double_peptide_final_data)
        
    if any(pep_calibrants_all['charge']==3):

        calibrants = triple_data
        accessions = deimos.get_accessions(peptide_data_file)
        data = deimos.load(peptide_data_file,
                           accession={'retention_time': 'MS:1000016',
                                      'drift_time': 'MS:1002476'})
        ms1 = data['ms1']
        ms1= ms1.drop(['function', 'process', 'scan'], axis=1)
        
        #nominal threshold
        ms1 = deimos.threshold(ms1, threshold=128)
        
        #sum over retention time
        ms1_sum = deimos.collapse(ms1, keep=['mz', 'drift_time'])
        ms1_sum= ms1_sum.sort_values(by='mz', ascending=True)
        
        calibrant_name = np.array(calibrants['calibrant'])
        mass = np.array(calibrants['m/z'])
        lit_ccs = np.array(calibrants['ccs'])
        
        extracted_drift_time = np.empty(shape= [0,1])
        
        def Gauss(x, amplitude, mean, stdev):
            y = amplitude * np.exp(-0.5*((x-mean)/stdev)**2)
            return y
        
        for i in range(0,len(mass)):
            ms_slice = deimos.slice(ms1_sum, by='mz', low=[mass[i]-0.05], high= [mass[i]+0.05])
            slice_peaks = deimos.peakpick.local_maxima(ms_slice, dims=['mz', 'drift_time'], bins=[37, 1])
            slice_peaks = slice_peaks.sort_values(by='drift_time', ascending=True)
            dt = np.array(slice_peaks['drift_time'])
            intensity = np.array(slice_peaks['intensity']) 
            parameters, covariance = curve_fit(Gauss, dt, intensity, maxfev=1000000, p0 = (max(intensity), dt[np.argmax(intensity)], 0.5))
            Gauss_intensity = round(parameters[0],0)
            Gauss_dt = round(parameters[1], 2)
            Gauss_stdev = round(parameters[2],2)
            extracted_drift_time = np.append(extracted_drift_time, Gauss_dt)
        
        # Calculating CCS Calibration Parameters
        
        corrected_dt = []
        corrected_ccs = []
        reduced_mass = []
        for i in range(0, len(mass)):
            reduced_mass_val = (mass[i] * nitrogen_mass)/(mass[i] + nitrogen_mass)
            reduced_mass = np.append(reduced_mass, reduced_mass_val)
            corrected_dt_val = extracted_drift_time[i] - (((math.sqrt(mass[i]))*edc)/1000)
            corrected_ccs_val = lit_ccs[i] * (math.sqrt(reduced_mass_val))/3
            corrected_dt = np.append(corrected_dt, corrected_dt_val)
            corrected_ccs = np.append(corrected_ccs, corrected_ccs_val)
        
        def cal_curve(corrected_dt, A, t0, B):
            y = A * (corrected_dt + t0)**B
            return y
        
        param, cov = curve_fit(cal_curve, corrected_dt, corrected_ccs, maxfev=1000000, p0=(500., 0.001, 0.5))
        
        triple_pep_A = param[0]
        triple_pep_t0 = param[1]
        triple_pep_B = param[2]
        
        # CALCULATES CCS VALUES AND RESIDUALS
        
        cal_ccs = []
        for i in range (0, len(corrected_dt)):
            cal_ccs_val = (triple_pep_A * math.sqrt(1/reduced_mass[i]) * (corrected_dt[i] + triple_pep_t0)**triple_pep_B)*3
            cal_ccs = np.append(cal_ccs, round(cal_ccs_val,2))
        
        residual = []
        for i in range(len(lit_ccs)):
            new_residual = round(((cal_ccs[i]-lit_ccs[i])/lit_ccs[i])*100,3)
            residual = np.append(residual, new_residual)
            
        # FINAL DATA OUTPUT
        
        triple_peptide_final_data = pd.DataFrame({'calibrant': calibrant_name, 'mz': mass, 'Drift Time': extracted_drift_time, 'Lit CCS': lit_ccs, 'Cal CCS': cal_ccs, 'Residuals': residual})
        triple_peptide_parameters = pd.DataFrame({'parameter': cal_parameter_names, 'values': [triple_pep_A, triple_pep_t0, triple_pep_B]})

experimental_data_exists = "exp_data_name" in locals()  

if experimental_data_exists: 
    exp_data = pd.read_excel(exp_data_name)
    compound = np.array(exp_data["Compound"])
    ion_mass = np.array(exp_data["m/z"])
    dt = np.array(exp_data["Drift time"])
    rt = np.array(exp_data["Retention time"])
    charge = np.array(exp_data["Charge"])
    
    cal_ccs = []
    corrected_dt = []
    biomolecule = []
    lipid_distances = []
    metabolite_distances = []
    single_pep_distances = []
    double_pep_distances = []
    triple_pep_distances = []
    lipid_ccs_values = []
    metabolite_ccs_values = []
    single_pep_ccs_values = []
    double_pep_ccs_values = []
    triple_pep_ccs_values = []
    biomolecules_present = [lipids_present, metabolites_present, peptides_present]
    
    for i in range(0, len(ion_mass)):
        reduced_mass = (ion_mass[i] * nitrogen_mass)/(ion_mass[i] + nitrogen_mass)
        corrected_dt = dt[i] - (((math.sqrt(ion_mass[i]))*edc)/1000)
        
        lipid_dist = 100000000
        metab_dist = 100000000
        single_peptide_dist = 100000000
        
        if lipids_present == 1:
            lipid_cal_ccs_val = round((lipid_A * math.sqrt(1/reduced_mass) * (corrected_dt + lipid_t0)**lipid_B), 3)
            lipid_line_CCS = ((ion_mass[i])**0.4788) * 11.837
            lipid_dist = np.linalg.norm((lipid_cal_ccs_val - lipid_line_CCS))
        if metabolites_present == 1:
            metab_cal_ccs_val = round((metab_A * math.sqrt(1/reduced_mass) * (corrected_dt + metab_t0)**metab_B), 3)
            metabolite_line_CCS = ((ion_mass[i])**0.3551) * 22.344
            metab_dist = np.linalg.norm((metab_cal_ccs_val - metabolite_line_CCS))
        if peptides_present == 1:
            if single_pep_A != 0:
                single_pep_cal_ccs_val = round((single_pep_A * math.sqrt(1/reduced_mass) * (corrected_dt + single_pep_t0)**single_pep_B), 3)
                single_peptide_line_ccs = ((ion_mass[i])**0.528) * 8.0136
                single_peptide_dist = np.linalg.norm((single_pep_cal_ccs_val - single_peptide_line_ccs))
            if double_pep_A != 0:
                double_pep_cal_ccs_val = round(((double_pep_A * math.sqrt(1/reduced_mass) * (corrected_dt + double_pep_t0)**double_pep_B)*2), 3)
            if triple_pep_A != 0:    
                triple_pep_cal_ccs_val = round(((triple_pep_A * math.sqrt(1/reduced_mass) * (corrected_dt + triple_pep_t0)**triple_pep_B)*3), 3)
        
        distances = [lipid_dist, metab_dist, single_peptide_dist]
        
        if all(biomolecules_present):
            if charge[i] == 2:
                cal_ccs = np.append(cal_ccs, double_pep_cal_ccs_val)
                assignment = "doubly charged peptide"
            elif charge[i] == 3:
                cal_ccs = np.append(cal_ccs, triple_pep_cal_ccs_val)
                assignment = "triply charged peptide"
            elif charge[i] == 1:
                if lipid_dist == min(distances):
                    if ion_mass[i] > 300: 
                        cal_ccs = np.append(cal_ccs, lipid_cal_ccs_val)
                        assignment = "lipid"
                    else:
                        cal_ccs = np.append(cal_ccs, metab_cal_ccs_val)
                        assignment = "small molecule"
                elif metab_dist == min(distances):
                    if ion_mass[i] < 550:
                        cal_ccs = np.append(cal_ccs, metab_cal_ccs_val)
                        assignment = "small molecule"
                    else:
                        cal_ccs = np.append(cal_ccs, single_pep_cal_ccs_val)
                        assignment = "singly charged peptide"
                elif single_peptide_dist == min(distances):
                    if ion_mass[i] > 300:
                        cal_ccs = np.append(cal_ccs, single_pep_cal_ccs_val)
                        assignment = "singly charged peptide"
                    else:
                       cal_ccs = np.append(cal_ccs, metab_cal_ccs_val)
                       assignment = "small molecule" 
                  
        #if only one biomolecule is present    
        elif lipids_present == 1 and (metabolites_present == 0 and peptides_present == 0):
            cal_ccs = np.append(cal_ccs, lipid_cal_ccs_val)
            assignment = "lipid"
                
        elif metabolites_present == 1 and (lipids_present == 0 and peptides_present == 0):
            cal_ccs = np.append(cal_ccs, metab_cal_ccs_val)
            assignment = "small molecule" 
            
        elif peptides_present == 1 and (lipids_present == 0 and metabolites_present == 0):
            if charge[i] == 1:
                cal_ccs = np.append(cal_ccs, single_pep_cal_ccs_val)
                assignment = "singly charged peptide"
            elif charge[i] == 2:
                cal_ccs = np.append(cal_ccs, double_pep_cal_ccs_val)
                assignment = "doubly charged peptide"
            elif charge[i] == 3:
                cal_ccs = np.append(cal_ccs, triple_pep_cal_ccs_val)
                assignment = "triply charged peptide"
            
        #if two biomolecules are present
        
        elif lipids_present == 1 and (metabolites_present == 1 and peptides_present == 0):
            if lipid_dist == min(distances):
                if ion_mass[i] > 300: 
                    cal_ccs = np.append(cal_ccs, lipid_cal_ccs_val)
                    assignment = "lipid"
                else:
                    cal_ccs = np.append(cal_ccs, metab_cal_ccs_val)
                    assignment = "small molecule"
            elif metab_dist == min(distances):
                if ion_mass[i] < 550:
                    cal_ccs = np.append(cal_ccs, metab_cal_ccs_val)
                    assignment = "small molecule"
                else:
                    cal_ccs = np.append(cal_ccs, lipid_cal_ccs_val)
                    assignment = "lipid"
            
        elif lipids_present == 1 and (metabolites_present == 0 and peptides_present == 1):
            if charge[i] == 2:
                cal_ccs = np.append(cal_ccs, double_pep_cal_ccs_val)
                assignment = "doubly charged peptide"
            elif charge[i] == 3:
                cal_ccs = np.append(cal_ccs, triple_pep_cal_ccs_val)
                assignment = "triply charged peptide" 
            elif charge[i] == 1:
                if lipid_dist == min(distances):
                    cal_ccs = np.append(cal_ccs, lipid_cal_ccs_val)
                    assignment = "lipid"
                elif single_peptide_dist == min(distances):
                    cal_ccs = np.append(cal_ccs, single_pep_cal_ccs_val)
                    assignment = "singly charged peptide" 
            
        elif lipids_present == 0 and (metabolites_present == 1 and peptides_present == 1):
            if charge[i] == 2:
                cal_ccs = np.append(cal_ccs, double_pep_cal_ccs_val)
                assignment = "doubly charged peptide"
            elif charge[i] == 3:
                cal_ccs = np.append(cal_ccs, triple_pep_cal_ccs_val)
                assignment = "triply charged peptide" 
            elif charge[i] == 1:
                if metab_dist == min(distances):
                    cal_ccs = np.append(cal_ccs, metab_cal_ccs_val)
                    assignment = "small molecule"
                elif single_peptide_dist == min(distances):
                    if ion_mass[i] > 300:
                        cal_ccs = np.append(cal_ccs, single_pep_cal_ccs_val)
                        assignment = "singly charged peptide"
                    else:
                       cal_ccs = np.append(cal_ccs, metab_cal_ccs_val)
                       assignment = "small molecule" 

        biomolecule = np.append(biomolecule, assignment)
        
        if lipids_present == 1 and (metabolites_present == 0 and peptides_present == 0):
            single_biomolecule = 1
        elif metabolites_present == 1 and (lipids_present == 0 and peptides_present == 0):
            single_biomolecule = 1
        elif peptides_present == 1 and (lipids_present == 0 and metabolites_present == 0):
            single_biomolecule = 1
        else:
            single_biomolecule = 0
            
# Calibration Effect        
        if single_biomolecule == 0:
            ccs_df = pd.DataFrame()
            if lipids_present == 1:
                lipid_ccs_values = np.append(lipid_ccs_values, lipid_cal_ccs_val)
                lipid_chosen_ccs = abs(lipid_ccs_values-cal_ccs)/cal_ccs*100
                ccs_df['Lipid CCS'] = lipid_chosen_ccs
            if metabolites_present == 1:
                metabolite_ccs_values = np.append(metabolite_ccs_values, metab_cal_ccs_val)
                metabolite_chosen_ccs = abs(metabolite_ccs_values-cal_ccs)/cal_ccs*100
                ccs_df['Metab CCS'] = metabolite_chosen_ccs
            if single_pep_A != 0:
                single_pep_ccs_values = np.append(single_pep_ccs_values, single_pep_cal_ccs_val)
                single_pep_chosen_ccs = abs(single_pep_ccs_values-cal_ccs)/cal_ccs*100
                ccs_df['P1 CCS'] = single_pep_chosen_ccs
    
            ccs_df['lowest'], ccs_df['second_lowest'] = np.sort(ccs_df, axis=1)[:, :2].T
    
            ccs_diff = np.array(ccs_df['second_lowest'])
            
            calibration_effect = []   
            for j in range(0,len(cal_ccs)):
                if charge[j] == 1:
                    if ccs_diff[j] <= 1:
                        calibration_effect = np.append(calibration_effect, 1)
                    elif ccs_diff[j] <= 3 and ccs_diff[j] > 1:
                        calibration_effect = np.append(calibration_effect, 2)
                    elif ccs_diff[j] <= 6 and ccs_diff[j] > 3:
                        calibration_effect = np.append(calibration_effect, 3)
                    else:
                        calibration_effect = np.append(calibration_effect, 4)
                else:
                    calibration_effect = np.append(calibration_effect, 0)
            
    if single_biomolecule == 1:
        final_data = pd.DataFrame({'Compound': compound, 'm/z': ion_mass, 'Charge': charge, 'Retention Time': rt, 'Drift Time': dt, 'Cal CCS': cal_ccs, 'Class Assignment': biomolecule})
    else :
        final_data = pd.DataFrame({'Compound': compound, 'm/z': ion_mass, 'Charge': charge, 'Retention Time': rt, 'Drift Time': dt, 'Cal CCS': cal_ccs,'Class Assignment': biomolecule, 'Calibration Effect': calibration_effect})

# Data Output    

writer = pd.ExcelWriter('Output/CalibrationOutput.xlsx', engine='xlsxwriter') 
if lipid_A != 0:         
    lipid_parameters.to_excel(writer, sheet_name='Lipid_Calibration', startrow = 0, startcol = 0)
    lipid_final_data.to_excel(writer, sheet_name='Lipid_Calibration', startrow = 6, startcol = 0)
if metab_A != 0:
    metabolite_parameters.to_excel(writer, sheet_name='SmallMolecule_Calibration', startrow = 0, startcol = 0)
    metabolite_final_data.to_excel(writer, sheet_name='SmallMolecule_Calibration', startrow = 6, startcol = 0)
if single_pep_A != 0:
    single_peptide_parameters.to_excel(writer, sheet_name='Peptide_z=1_Calibration', startrow = 0, startcol = 0)
    single_peptide_final_data.to_excel(writer, sheet_name='Peptide_z=1_Calibration', startrow = 6, startcol = 0)
if double_pep_A != 0:
    double_peptide_parameters.to_excel(writer, sheet_name='Peptide_z=2_Calibration', startrow = 0, startcol = 0)
    double_peptide_final_data.to_excel(writer, sheet_name='Peptide_z=2_Calibration', startrow = 6, startcol = 0)
if triple_pep_A != 0:
    triple_peptide_parameters.to_excel(writer, sheet_name='Peptide_z=3_Calibration', startrow = 0, startcol = 0)
    triple_peptide_final_data.to_excel(writer, sheet_name='Peptide_z=3_Calibration', startrow = 6, startcol = 0)
if experimental_data_exists:
    final_data.to_excel(writer, sheet_name='Experimental_Output', startrow = 0, startcol = 0)
writer.close()

ws = Tk()
ws.title('Processing Complete!')
ws.geometry('400x100') 

def nextPage():
    ws.destroy()

Label(ws, text='Check Output Folder for Results', foreground='green').place(relx = 0.5, rely = 0.5, anchor= CENTER)
 
next_btn = Button(ws, text ='OK', command = nextPage).place(relx = 0.5, rely = 1.0, anchor= S)

ws.mainloop()
   
#GUI CODE: Developed from PythonGuides
import pandas as pd
import tkinter as tk
from tkinter import filedialog
import os

#%% 
application_window = tk.Tk() 
    
fTyp = [("fichier de données (*.csv)", "*.csv")]
#Ask the user to choose a file to convert
TB_file_name = filedialog.askopenfilenames(parent=application_window,
                                    initialdir=os.getcwd(),
                                    title="Please select one or more file(s) csv file containing the Draw tab:",
                                    filetypes=fTyp)


application_window.destroy()

application_window = tk.Tk()
destinationPath = filedialog.askdirectory(parent=application_window,
                                 initialdir=os.getcwd(),
                                 title="Please select a folder to put your translated dataFrame(s):")
application_window.destroy()

#%%
#We'll start the drawing 15 min after the beginning 
for n in range(len(TB_file_name)):
    TB_Draw_Tab = pd.read_csv(TB_file_name[n])   
    Sim_Draw_Tab = pd.DataFrame( columns = ['Start(h)', 'Duration(min)', 'Draw Rate(L/min)'])
    
    #Add the first line
    new_row = pd.DataFrame({"Start(h)": 15/60,\
                         "Duration(min)": TB_Draw_Tab['Duration Modified (s)'][0]/60,\
                         "Draw Rate(L/min)": TB_Draw_Tab["Modified flow rate (L/min)"][0]}, index = [0])
    Sim_Draw_Tab = pd.concat([Sim_Draw_Tab, new_row], ignore_index= True)
    I = len(TB_Draw_Tab.index)
    #Add the other rows
    for i in range(1,I):
        new_row = pd.DataFrame({"Start(h)": (TB_Draw_Tab['Wait Time (s)'][i-1] + TB_Draw_Tab['Duration Modified (s)'][i-1])/3600,\
                              "Duration(min)": TB_Draw_Tab['Duration Modified (s)'][0]/60,\
                              "Draw Rate(L/min)": TB_Draw_Tab["Modified flow rate (L/min)"][i]}, index = [0])
        Sim_Draw_Tab = pd.concat([Sim_Draw_Tab, new_row], ignore_index= True)
        
#%%Now saving the dataFrame in a file
#destination folder
    SimFileName = os.path.basename(TB_file_name[n])[:-4]
    
    
    SimFileName = destinationPath + '/' +  SimFileName + '_Sim.csv'
    
    Sim_Draw_Tab.to_csv(SimFileName)























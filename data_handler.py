import pandas as pd
import numpy as np

# data_dict={"param": [1, 2, 3, [4,5,6]],
#                         "model"  :   "Shu",
#                         "radius" :   [[1,2,3], [1,2,3]],
#                         "density":   [],
#                         "v_r"    :   [],
#                         "v_theta":   [],
#                         "v_phi"  :   []}

# mid_term_marks_df = pd.DataFrame(data_dict)
# mid_term_marks_df.to_csv("midterm.csv")
# print(mid_term_marks_df)
df = pd.read_csv("data.csv")
print(df)

class data():
    def __init__(self):
        self.df = pd.read_csv("data.csv")
    
    def save(self, param=[], model="", radius=[], density=[], v_r=[], v_theta=[], v_phi=[]):
        """
        Saving data into csv. list -> csv

        args:
            param: the parameter used to calculate data
                cs, t_acc, omega, time 
        """        
        ifSave = self.isSaved(param, model)

        if ifSave == False:
            print("Saving data...")
            data = [param, model, radius, density, v_r, v_theta, v_phi]
            self.df.loc[len(self.df)] = data            
            self.df.to_csv('data.csv')
    
    def load(self, param=[], model=""):
        tf_param = list(self.df['param'].isin([str(param)]))
        tf_model = list(self.df['model'].isin([str(model)]))
        print(tf_param, tf_model)
        for i in range(len(tf_param)):
            if tf_param[i] == True and tf_model[i] == True:
                radius  = self.df.iloc[i]['radius']
                density = self.df.iloc[i]['radius']
                v_r     = self.df.iloc[i]['v_r']
                v_theta = self.df.iloc[i]['v_theta']
                v_phi   = self.df.iloc[i]['v_phi']

        return radius, density, v_r, v_theta, v_phi
    
    def isSaved(self, param=[], model=""):
        """
        Check if the combination of parameter and model have been save enough.

        args:
            param: the parameter used to calculate data
                [cs, t_acc, omega, time]
            model: which version of calculation
                "Shu" or "TSC"

        return:
            ifSaved: determine if the combination have save before
                True: yes, the data have been saved before
                False: no, the data need to be save
        """
        ifSave = False
        tf_param = list(self.df['param'].isin([str(param)]))
        tf_model = list(self.df['model'].isin([str(model)]))
        print(tf_param, tf_model)
        for i in range(len(tf_param)):
            if tf_param[i] == True and tf_model[i] == True:
                ifSave = True
                break
        
        if ifSave == True:
            print("Data already exist!")
        
        return ifSave
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
            print("Data already exist!")

        if ifSave == True:
            print("Saving data...")
            data = [param, model, radius, density, v_r, v_theta, v_phi]
            self.df.loc[len(self.df)] = data            
            self.df.to_csv('data.csv')
    
    def isSaved(self, param=[], model=""):
        ifSave = True
        tf_param = list(self.df['param'].isin([str(param)]))
        tf_model = list(self.df['model'].isin([str(model)]))
        print(tf_param, tf_model)
        for i in range(len(tf_param)):
            if tf_param[i] == True and tf_model[i] == True:
                ifSave = False
                break
        
        return ifSave
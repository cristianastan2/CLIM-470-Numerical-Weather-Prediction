import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
from datetime import datetime

data = pd.read_csv('Book2.csv')
data.columns = ['Day', 'tmin(Forecast)', 'tmin(Verification)']

for col in ['tmin(Forecast)', 'tmin(Verification)']:
    data[col] = pd.to_numeric(data[col])
    
data =  data.set_index('Day')
plt.figure(figsize=(12,6))
plt.xlabel('Day', fontsize=15)
plt.ylabel('Temperature (F)', fontsize=12)
plt.plot(data.index, data['tmin(Forecast)'])
plt.plot(data.index, data['tmin(Verification)'])

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px


dirPath = "../data"

if __name__ == "__main__":
    df = pd.read_csv(f"{dirPath}/check_convection.txt", names=['field', 'x', 'y', 'val', 'val2'], delimiter=" ", index_col=False) 
    df.replace("X", "Hx", inplace=True)
    df.replace("Y", "Hy", inplace=True)

    df.loc[df['field']=='Hx', 'x'] += 0.001
    df.loc[df['field']=='Hx', 'y'] += 0.001
    df.loc[df['field']=='Hy', 'x'] += 0.001
    df.loc[df['field']=='Hy', 'y'] += 0.001
    # df = df.apply(pd.to_numeric, args=('coerce',))

    print(df[df["field"] == "Hx"].head(5))
    df_u = df[df['field'] == 'u']
    df_v = df[df['field'] == 'v']
    df_Hx = df[df['field'] == 'X']
    df_Hy = df[df['field'] == 'Y']


    fig = px.scatter(df, x='x', y='y', color='field', symbol='field', hover_data=["val", "val2"])
    fig.add_hline(y=3.)
    fig.add_hline(y=2.)
    fig.add_vline(x=3.)
    fig.add_vline(x=8.)

    fig.show()


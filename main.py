import csv
import numpy as np
import pandas as pd
import os


temp_output = []
co2_column_output = []

for mcd_file in os.listdir("./MCD"):
    mcd_df = pd.read_csv(f"./MCD/{mcd_file}", sep=";")

    result = {"sol": int(mcd_file.split(".")[0].split('_')[-2])}
    for simu_file in os.listdir("./simulation_results/sim1"):
        simu_df = pd.read_csv(f"./simulation_results/sim1/{simu_file}")

        simu_temperature = simu_df["n.temp"]
        simu_co2_colummn = simu_df["n.co2_column"]

        mcd_temperature = mcd_df["temperture"]
        mcd_co2_column = mcd_df["extvar_66"]

        temp_cov = np.corrcoef(simu_temperature, mcd_temperature)[0][1]
        co2_column_cov = np.corrcoef(simu_co2_colummn, mcd_co2_column)[0][1]

        result[simu_file] = {"temp": temp_cov, "co2_column": co2_column_cov}

    temp_dict = dict()
    co2_column_dict = dict()

    for k, v in result.items():
        if k == "sol":
            temp_dict[k] = v
            co2_column_dict[k] = v
        else:
            temp_dict[k] = v["temp"]
            co2_column_dict[k] = v["co2_column"]

    temp_output.append(temp_dict)
    co2_column_output.append(co2_column_dict)

temp_cov_df = pd.DataFrame(temp_output)
co2_column_cov_df = pd.DataFrame(co2_column_output)

temp_cov_df.sort_values(by=['sol']).to_csv("./results/temp_cov.csv", index=False)
co2_column_cov_df.sort_values(by=['sol']).to_csv("./results/co2_column_cov.csv", index=False)


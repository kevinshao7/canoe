import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata
from tqdm import tqdm
import os
from os.path import isfile, join
import pathlib


# Path to your single combined .nc file
case = "uranus2m"
mypath = str(pathlib.Path(__file__).parent.resolve())+"/"
# Variables to process
variables_to_process2 = ["vel1","press"]
variables_to_process3 = ["temp", "theta", "vlat", "vlon"]

num_time_steps=300
timearr=[]

def process_timestep(t):
    filename ="{}.out2.{:05d}.nc".format(case,t)
    with Dataset(filename, mode="r") as nc:
        # x1 = nc.variables["x1"][:]
        # x2 = nc.variables["x2"][:]
        # x3 = nc.variables["x3"][:]
        timearr.append(nc.variables["time"][:])
        data_arrays2 = {
            var: nc.variables[var][:, :, :] for var in variables_to_process2
        }
    filename = "{}.out3.{:05d}.nc".format(case,t)
    with Dataset(filename, mode="r") as nc:
        # x1 = nc.variables["x1"][:]
        # x2 = nc.variables["x2"][:]
        # x3 = nc.variables["x3"][:]
        # time_value = nc.variables["time"][:]

        data_arrays3 = {

            var: nc.variables[var][:, :, :] for var in variables_to_process3
        }

    return data_arrays2,data_arrays3


# Process each time step and store results in a list
results2 = []
results3 = []
for t in tqdm(range(1,num_time_steps), desc="Processing time steps"):
    arr2,arr3 = process_timestep(t)
    results2.append(arr2)
    results3.append(arr3)

firstfile= "{}.out2.00001.nc".format(case)
with Dataset(firstfile, "r") as nc:
    x1 = nc.variables["x1"][:]
    x2 = nc.variables["x2"][:]
    x3 = nc.variables["x3"][:]

# Save the processed data to a new NetCDF file
output_filepath = mypath+case+"_main.nc"
with Dataset(output_filepath, mode="w") as new_nc:
    # Create the dimensions
    new_nc.createDimension("time", None)  # Unlimited dimension (usually time)
    new_nc.createDimension("x1", len(x1))
    new_nc.createDimension("x2", len(x2))
    new_nc.createDimension("x3", len(x3))

    # Create and assign the time variable
    time_var = new_nc.createVariable("time", np.float64, ("time",))
    time_var[:] = timearr

    # Create and assign the x1 variable
    new_nc.createVariable("x1", np.float32, ("x1",))
    new_nc.variables["x1"][:] = x1

    new_nc.createVariable("x2", np.float32, ("x2",))
    new_nc.variables["x2"][:] = x2

    new_nc.createVariable("x3", np.float32, ("x3",))
    new_nc.variables["x3"][:] = x3
    # Create and assign the processed data
    for var in variables_to_process2:
        new_var = new_nc.createVariable(var, np.float32, ("time", "x1", "x2", "x3"))
        for t in range(num_time_steps):
            new_var[t, :, :, :] = results2[t-1][var]
        # Create and assign the latitude and longitude variables
    for var in variables_to_process3:
        new_var = new_nc.createVariable(var, np.float32, ("time", "x1", "x2", "x3"))
        for t in range(num_time_steps):
            new_var[t, :, :, :] = results3[t-1][var]
        # Create and assign the latitude and longitude variables

print(f"Processing completed. Results saved to {output_filepath}.")

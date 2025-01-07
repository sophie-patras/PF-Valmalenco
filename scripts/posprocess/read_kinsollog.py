import re
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('/home/patras/PF-Valmalenco/scripts/config.mplstyle')
#path_fig = '/mnt/c/Users/User/Documents/POLIMI/0_TESI/8-Figures/'   #local
path_fig = '/mnt/c/Users/Sophie/Documents/4-Figures/'   #distant
path_fig = '/mnt/c/Users/User/Documents/POLIMI/0_TESI/8-Figures/'   #distant
###

def parse_log_file(log_file_path):
    # Array to hold the datasets

    # Regular expressions for the key data
    time_step_pattern = re.compile(r"KINSOL starting step for time (\d+\.\d+)")
    nonlin_its_pattern = re.compile(r"Nonlin. Its.:.*?(\d+)")
    lin_its_pattern = re.compile(r"Lin. Its.:.*?(\d+)")
    func_evals_pattern = re.compile(r"Func. Evals.:.*?(\d+)")
    pc_evals_pattern = re.compile(r"PC Evals.:.*?(\d+)")
    pc_solves_pattern = re.compile(r"PC Solves:.*?(\d+)")
    lin_conv_fails_pattern = re.compile(r"Lin. Conv. Fails:.*?(\d+)")
    beta_cond_fails_pattern = re.compile(r"Beta Cond. Fails:.*?(\d+)")
    backtracks_pattern = re.compile(r"Backtracks:.*?(\d+)")

    with open(log_file_path, 'r') as file:
        content = file.read()

    # Find all time steps and other attributes in the log file
    time_steps = time_step_pattern.findall(content)
    tot_ite = len(time_steps)
    #print('nb of iteration',tot_ite)

    dataset = np.empty((tot_ite,9))
    dataset[:,0] = time_steps
    dataset[:,1] = nonlin_its_pattern.findall(content)
    dataset[:,2] = lin_its_pattern.findall(content)
    dataset[:,3] = func_evals_pattern.findall(content)
    dataset[:,4] = pc_evals_pattern.findall(content)
    dataset[:,5] = pc_solves_pattern.findall(content)
    dataset[:,6] = lin_conv_fails_pattern.findall(content)
    dataset[:,7] = beta_cond_fails_pattern.findall(content)
    dataset[:,8] = backtracks_pattern.findall(content)

    return dataset

# Example usage:
#log_file_path = "path_to_your_log_file.txt"

path = '/home/patras/PF-Valmalenco/outputs/'
foldername = 'CLM_V2'
runname = 'CLM_V2'
log_file_path = f'{path}{foldername}/{runname}.out.kinsol.log'

parsed_data = parse_log_file(log_file_path)

# To print the parsed data
#for data in parsed_data:
print(parsed_data)

# print nb of nni for each time step
fig = plt.figure(1)
ax = fig.add_subplot(111)

ax.plot(parsed_data[:,0],parsed_data[:,1],marker='+',linestyle='',label='nli')
ax.plot(parsed_data[:,0],parsed_data[:,2],marker='+',linestyle='',label='li')
ax.grid(True)
ax.legend()
ax.set_xlabel('time step [h]')
ax.set_ylabel('number of iterations')

plt.savefig(f'{path_fig}{foldername}.nni.kinsol.png')
import re
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('/home/patras/PF-Valmalenco/scripts/config.mplstyle')
#path_fig = '/mnt/c/Users/User/Documents/POLIMI/0_TESI/8-Figures/'   #local
path_fig = '/mnt/c/Users/Sophie/Documents/4-Figures/'   #distant
#path_fig = '/mnt/c/Users/User/Documents/POLIMI/0_TESI/8-Figures/'   #distant
###

#########################################################################
"""
# kinsol type by time step

KINSOL starting step for time <simulated_time:time>
scsteptol used:  <Solver.Nonlinear.StepTol, es. 1e-30>  
fnormtol  used:  <Solver.Nonlinear.ResidualTol, es. 1e-06>
% for every nonlinear iteration (nni), this is listing the residual (fnorm) and the
% number of function evaluations (nfe) - ie. k iteration of the suite ? - that where
% completed within this iteration.
KINSolInit nni=    0  fnorm=        0.02623070719595643  nfe=     1 
KINSol nni=    1 fnorm=        0.02619648992107566 nfe=     3
...
KINSol nni=    <x<Nonlinear.MaxIter> fnorm=        0.02046724352021377 nfe=   166
KINSol return value <1,3,4>
---KINSOL_<SUCCESS,STEP_LT_STPTOL,LNSRCH_NONCONV>
% success, fnorm stagnated with difference below StepTol, exceed nb MaxIter
% In this 2 last case, PF will try solving with half step for <Solver.MaxConvergenceFailures> loops maximum
\n
-------------------------------------------------- 
                    Iteration             Total
Nonlin. Its.:              24                24
Lin. Its.:                 26                26
Func. Evals.:            1591              1591
PC Evals.:                 24                24
PC Solves:                 50                50
Lin. Conv. Fails:           0                 0
Beta Cond. Fails:           1                 1
Backtracks:                96                96
-------------------------------------------------- 
\n
"""
#########################################################################

def parse_log_file(log_file_path):
    # Array to hold the datasets

    # Regular expressions for the key data
    time_step_pattern = re.compile(r"KINSOL starting step for time (\d+\.\d+)")
    successvalue_pattern = re.compile(r'KINSol return value *?(\d+)')
    success_pattern = re.compile(r"---KINSOL_\w+")
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

    dataset = np.empty((tot_ite,10))
    dataset[:,0] = time_steps
    dataset[:,1] = nonlin_its_pattern.findall(content)
    dataset[:,2] = lin_its_pattern.findall(content)
    dataset[:,3] = func_evals_pattern.findall(content)
    dataset[:,4] = pc_evals_pattern.findall(content)
    dataset[:,5] = pc_solves_pattern.findall(content)
    dataset[:,6] = lin_conv_fails_pattern.findall(content)
    dataset[:,7] = beta_cond_fails_pattern.findall(content)
    dataset[:,8] = backtracks_pattern.findall(content)

    idx_NOKIN = np.where(dataset[:,1]==0)
    #print(idx_NOKIN)
    dataset[idx_NOKIN, 9] = 2.0
    idx_else = [i for i in range(tot_ite) if i not in idx_NOKIN[0]]
    dataset[idx_else,9] = successvalue_pattern.findall(content)

    #if len(successvalue_pattern.findall(content)) == tot_ite:
    #    dataset[:,9] = successvalue_pattern.findall(content)

    #successlist = success_pattern.findall(content)
    #print(len(successlist ))
    print("Failure at sequence(s):") #, len(dataset[:,9])) #successlist)
    #print(np.where(successlist != '---KINSOL_SUCCESS'))
    print(np.where(np.array(dataset[:,9]) != 1.0)) # ie. idx_NOKIN

    return dataset

# Example usage:
#log_file_path = "path_to_your_log_file.txt"

path = '/home/patras/PF-Valmalenco/outputs/'
foldername = 'CLM_V52'
runname = 'CLM_V5'

## MX
#path = '/home/patras/PF-Test/Maxwell2013/outputs/'
#foldername = 'MX.c1s1y3_v6'
#runname = 'MX.c1s1'

## DS
#path = '/home/patras/PF-Test/DumbSquare/outputs/'
#foldername = 'DSc100z10s0.MX.DP23.IN0_v48' #'DS.c100s1_v28'
#runname = 'DSc100z10s0'

log_file_path = f'{path}{foldername}/{runname}.out.kinsol.log'

parsed_data = parse_log_file(log_file_path)

# To print the parsed data
#for data in parsed_data:
#print(parsed_data)

shape_data = parsed_data.shape
nt = shape_data[0]

## final time step to plot
tfin = -1
#tfin = 136

idx_SUCCESS = np.where(np.array(parsed_data[:tfin,9]) == 1.0)
#print(idx_SUCCESS)
idx_STEP_LT_STPTOL = np.where(np.array(parsed_data[:tfin,9]) == 3.0)
#print(idx_STEP_LT_STPTOL)
idx_LNSRCH_NONCONV = np.where(np.array(parsed_data[:tfin,9]) == 4.0)
# Use 2 when KINSOL is not evaluated
idx_NOKIN = np.where(np.array(parsed_data[:tfin,9]) == 2.0) # Perso >> NO FLAG FOR 2


# print nb of nni for each time step
fig = plt.figure(1,figsize=(5,3))
ax = fig.add_subplot(111)


if len(idx_SUCCESS[0]) < 2:
    # GLOBAL
    ax.plot(parsed_data[:tfin,0],parsed_data[:tfin,1],marker='+',linestyle='-',label='nli')
    ax.plot(parsed_data[:tfin,0],parsed_data[:tfin,2],marker='+',linestyle='-',label='li')
else:
    # SPECIFIC FLAGS
    ax.plot(parsed_data[idx_SUCCESS,0],parsed_data[idx_SUCCESS,1],marker='+',linestyle='-',color='b',label='nli-SUCCESS')
    ax.plot(parsed_data[idx_SUCCESS,0],parsed_data[idx_SUCCESS,2],marker='+',linestyle='-',color='g',label='li-SUCCESS')
    ax.plot(parsed_data[idx_STEP_LT_STPTOL,0],parsed_data[idx_STEP_LT_STPTOL,1],marker='s',color='b',linestyle='-',label='nli-STEP_LT_STPTOL')
    ax.plot(parsed_data[idx_STEP_LT_STPTOL,0],parsed_data[idx_STEP_LT_STPTOL,2],marker='s',color='g',linestyle='-',label='li-STEP_LT_STPTOL')
    ax.plot(parsed_data[idx_LNSRCH_NONCONV,0],parsed_data[idx_LNSRCH_NONCONV,1],marker='*',color='b', linestyle='-',label='nli-LNSRCH_NONCONV')
    ax.plot(parsed_data[idx_LNSRCH_NONCONV,0],parsed_data[idx_LNSRCH_NONCONV,2],marker='*',color='g', linestyle='-',label='li-LNSRCH_NONCONV')
    ax.plot(parsed_data[idx_NOKIN,0],parsed_data[idx_NOKIN,1],marker='x',color='b', linestyle='-',label='nli-NOKIN')
    ax.plot(parsed_data[idx_NOKIN,0],parsed_data[idx_NOKIN,2],marker='x',color='g', linestyle='-',label='li-NOKIN')


ax.grid(True)
#ax.legend()
#ax.set_ylim((0,100))
ax.set_xlabel('time [h]')
ax.set_ylabel('number of iterations')

plt.savefig(f'{path_fig}{foldername}.ni.kinsol.png')
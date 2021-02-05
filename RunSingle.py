'''
    RunSingle
    Start a simulation of the agent-based model by launching cpp-code 'swarmdyn'
    with the here defined parameters.
    Also visualizes the simulation by calling the script "AnimateRun.py".
    v0.1, 13.5.2020

    (C) 2020 Pascal Klamser

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import os
import time as pytime
import matplotlib
if __name__ == '__main__':
    matplotlib.use('TkAgg')
import SwarmDynByPy as swarmPy
import AnimateRun
import json
from pathlib import Path

def main():
    # Input Parameters
    #########################################
    dockerName = None # alternatively 'gcc_docker' or None see README.md for usage (Docker-alternative)

    trans_time = 100 # 200
    record_time = 40 # 1200

    dic = dict()
    dic = swarmPy.get_base_params(record_time, trans_time=trans_time)
    # changes from base-parameters
    dicC = dict()
    dicC['N'] = 400
    dicC['output_mode'] = 2
    # dic['output'] = 0.2
    dicC['beta'] = 0.4 # relaxation coefficient
    dicC['alg_strength'] = 9.9 # the order-disorder transition is at alg_strength=0.83
    # dicC['rep_strength'] = 2.26 # the order-disorder transition is at alg_strength=0.83
    f_name = ''
    keys = sorted(list(dicC.keys()))
    for k in keys:
        f_name += '{}{}'.format(k, dicC[k])
        dic[k] = dicC[k]
    print(f_name)


    # Generate and Run Command
    #########################################
    swarmPy.solve_parameter_dependencies(dic)
    f_para = Path.cwd() / 'parameters.json'
    json.dump(dic, f_para.open('w', encoding='utf-8'), indent=4, ensure_ascii=False)
    command = swarmPy.dic2swarmdyn_command(dic)
    print(command)
    t0 = pytime.time()
    dockerCall = ''
    if dockerName is not None:
        dockerCall = 'docker run --rm -v "$PWD":/usr/src/myapp -w /usr/src/myapp {} '.format(dockerName)
    os.system(dockerCall + 'make cl;')
    os.system(dockerCall + command)
    t1 = pytime.time()
    print(t1-t0)
    AnimateRun.main(None)
    return

if __name__ == '__main__':
    main()

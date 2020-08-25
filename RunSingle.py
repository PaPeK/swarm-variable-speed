'''
    RunSingle
    Start a simulation of the agent-based predator prey  model by launching cpp-code 'swarmdyn'
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
import SwarmDynByPy as swarmPy
import AnimateRun

def main():
    # Input Parameters
    #########################################
    dockerName = None # alternatively 'gcc_docker' or None see README.md for usage (Docker-alternative)

    pred_time = 200
    record_time = 60

    dic = dict()
    dic = swarmPy.get_base_params(pred_time, record_time)
    # changes from base-parameters
    dic['pava_sig'] = 0.1 # default 0.075
    dic['alg_strength'] = 3 # the order-disorder transition is at alg_strength=0.83

    # Generate and Run Command
    #########################################
    pavas = swarmPy.heteroPopulation(dic) # to change heterogeneity: modify dic['pava_sig']
    swarmPy.solve_parameter_dependencies(dic)
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

#!/usr/bin/env python3

import stim
import sys
import copy
import networkx as nx
import math
import os
import cairosvg
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pymatching
import random
matplotlib.get_backend()

class Qubit(object):
    def __init__(self, x, y, ind):
        self.x = x
        self.y = y
        self.ind = ind
        self.processed = False
        self.isancilla = False
        self.type = "NA"

    def getind(self):
        return self.ind

    def isprocessed(self):
        return self.processed
    
    def setprocessed(self):
        self.processed = True

    def setancilla(self, t):
        self.isancilla = True
        self.type = t

    def getancillatype(self):
        return self.type

    def gety(self):
        return self.y

    def getx(self):
        return self.x

    def getancillatype(self):
        return self.type

    def printme(self, pref=""):
        print(F'{pref} {self.x} {self.y} {self.ind} {self.processed}') 

class Vertex(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def addqbitinds(self, qbitinds):
        self.qubitinds = qbitinds

    def getqbitinds(self):
        return self.qubitinds 

    def getx(self):
        return self.x
    
    def gety(self):
        return self.y

    def addancillaqbitind(self, ind):
        self.ancillaqbitind = ind

    def getancillaind(self):
        return self.ancillaqbitind

class Plaquette(object):
    def __init__(self, locx, locy):
        self.qubitinds = []
        self.x = locx
        self.y = locy

    def addqbitinds(self, qbitinds):
        self.qubitinds = qbitinds

    def addancillaqbitind(self, ind):
        self.ancillaqbitind = ind

    def getancillaind(self):
        return self.ancillaqbitind

    def getqbitinds(self):
        return self.qubitinds

    def printme(self, pref=""):
        print(F'{pref} {self.x} {self.y}')

    def getx(self):
        return self.x

    def gety(self):
        return self.y

class Topology(object):
    def __init__(self, k, r, c):
        self.nplaquettes = k
        self.rows = r
        self.columns = c
        self.plaquettes = {}
        self.qinddict = {}
        self.vertexinddict = {}
        self.vertices = {}
        self.qubits = {}
        self.ancillaqubits = {}
        self.sqbitindmap = {}
        self.qind = 0
        self.initialized = False
        self.initializetopology()
        self.noise = -1
        self.circuit = stim.Circuit()
        
    def setNoise(self, noise):
        self.noise = noise

    def addNoise(self, qbitinds):
        if self.noise > 0:
            if len(qbitinds) > 1:
                self.circuit.append_operation("DEPOLARIZE2", qbitinds, self.noise)
            elif len(qbitinds) > 0:
                self.circuit.append_operation("DEPOLARIZE1", qbitinds, self.noise)

    def getNumQubits(self):
        return len(self.qubits)

    def getqbits(self, i, j):
        xs = [i - 1, i, i + 1, i]
        ys = [j, j - 1, j, j + 1]
        rqubitinds = []
        for xbar, ybar in zip(xs, ys):
            x = xbar%(2*self.rows)
            y = ybar%(2*self.columns)
            if "{}-{}".format(x, y) in self.qinddict:
                rqubitinds.append(self.qinddict["{}-{}".format(x, y)])
            else:
                qbit = Qubit(x, y, self.qind)
                self.qinddict["{}-{}".format(x,y)] = self.qind 
                self.qubits[self.qind] = qbit
                self.sqbitindmap[self.qind] = self.qind
                rqubitinds.append(self.qind)
                self.qind += 1
        return rqubitinds

    def getvertexqbits(self, xv, yv):
        xs = [xv-1, xv, xv+1, xv]
        ys = [yv, yv-1, yv, yv+1]
        rqubitinds = []
        for xbar, ybar in zip(xs, ys):
            x = xbar%(2*self.rows)
            y = ybar%(2*self.columns)
            for qind in self.qubits:
                qx = self.qubits[qind].getx()
                qy = self.qubits[qind].gety()
                if qx == x and qy == y:
                    rqubitinds.append(qind)
        return rqubitinds

    def initializevertices(self):
        vind = 0
        for i in range(self.nplaquettes):
            xi = self.plaquettes[i].getx()
            yi = self.plaquettes[i].gety()
            xs = [xi-1,xi-1,xi+1,xi+1]
            ys = [yi+1,yi-1,yi-1,yi+1]
            for xbar,ybar in zip(xs,ys):
                x = xbar%(2*self.rows)
                y = ybar%(2*self.columns)
                if "{}-{}".format(x, y) not in self.vertexinddict:
                    self.vertexinddict["{}-{}".format(x,y)] = vind
                    vertex = Vertex(x, y)
                    vertex.addqbitinds(self.getvertexqbits(x, y))
                    self.vertices[vind] = vertex
                    vind += 1

    def initializetopology(self):
        ind = 0
        x = 1
        y = 1
        for i in range(self.rows):
            for j in range(self.columns):
                self.plaquettes[ind] = Plaquette(x, y)
                self.plaquettes[ind].addqbitinds(self.getqbits(x, y))
                ind += 1
                if i % 2 == 0:
                    y += 2
                else:
                    y -= 2
            x += 2
            if i % 2 == 0:
                y = self.columns*2 - 1 
            else:
                y = 1 
        self.initializevertices()
        self.initialized = True

    def implementcircuit(self):
        if not self.initialized:
            print("ERROR topology has not been initialized")
            return
        allinds = []
        for qind in self.qubits:
            self.circuit.append_operation("QUBIT_COORDS", [self.sqbitindmap[qind]], [self.qubits[qind].gety(), self.qubits[qind].getx()] )
            allinds.append(self.sqbitindmap[qind])
        self.circuit.append_operation("R", allinds)
        self.circuit.append_operation("TICK")
        for i in range(self.nplaquettes-1):
            qbitinds = copy.deepcopy(self.plaquettes[i].getqbitinds())
            cqbitind = -1
            for qi in qbitinds:
                if self.qubits[qi].isprocessed() == False:
                    cqbitind = qi
                    break
            qbitinds.remove(cqbitind)
            self.qubits[cqbitind].setprocessed()
            self.circuit.append_operation("H", self.sqbitindmap[cqbitind])
            self.addNoise([self.sqbitindmap[cqbitind]])
            for qi in qbitinds:
                self.circuit.append_operation("CNOT", [self.sqbitindmap[cqbitind], self.sqbitindmap[qi]])
                self.addNoise([self.sqbitindmap[cqbitind], self.sqbitindmap[qi]])
                self.qubits[qi].setprocessed()

    def addancilla(self, check='p'):
        self.ancillaqubits = {}
        if check=='p':
            for i in self.plaquettes:
                xi = self.plaquettes[i].getx()
                yi = self.plaquettes[i].gety()
                qbit = Qubit(xi, yi, self.qind)
                qbit.setancilla("p")
                self.qinddict["{}-{}".format(xi,yi)] = self.qind 
                self.ancillaqubits[self.qind] = qbit
                self.sqbitindmap[self.qind] = self.qind
                self.plaquettes[i].addancillaqbitind(self.qind)
                self.qind += 1
        if check=='v':
            for vind in self.vertices: 
                xvi = self.vertices[vind].getx()
                yvi = self.vertices[vind].gety()
                qbit = Qubit(xvi, yvi, self.qind)
                qbit.setancilla("v")
                self.qinddict["{}-{}".format(xvi,yvi)] = self.qind 
                self.ancillaqubits[self.qind] = qbit
                self.sqbitindmap[self.qind] = self.qind
                self.vertices[vind].addancillaqbitind(self.qind)
                self.qind += 1
        for qind in self.ancillaqubits:
            self.circuit.append_operation("QUBIT_COORDS", [self.sqbitindmap[qind]], [self.ancillaqubits[qind].gety(), self.ancillaqubits[qind].getx()] )
        ancillainds = sorted(self.ancillaqubits.keys())
        self.circuit.append_operation("MR", ancillainds)
        self.circuit.append_operation("TICK")

    def printcircuit(self, visual='timeline'):
        if not self.initialized:
            print("ERROR topology has not been initialized")
            return
        if visual == 'None':
            print(self.circuit)
            return
        if visual == 'simple':
            print(self.circuit.diagram())
            return
        if visual == 'timeline':
            ckttype = 'timeline-svg'
        if visual == 'timeslice':
            ckttype = 'timeslice-svg'
        if os.path.exists('out.svg'):
            os.remove('out.svg')
        with open('out.svg', 'w') as fout:
            for l in str(self.circuit.diagram(ckttype)).split('\n'):
                fout.write(l)
        cairosvg.svg2png(url='out.svg', write_to='out.png')
        imageobj = plt.imread('out.png')
        plt.imshow(imageobj)
        plt.show()
        return

    def addmeasurement(self):
        if not self.initialized:
            print("ERROR topology has not been initialized")
            return
        self.circuit.append_operation("MR", self.qubits.keys())
        
    def takeameasurement(self, detector=False, num_shots=500):
        if not self.initialized:
            print("ERROR topology has not been initialized")
            return
        if not detector:
            samples = self.circuit.compile_sampler().sample(1)[0]
        else:
            samples = self.circuit.compile_detector_sampler().sample(shots=num_shots, separate_observables=True)
        return samples

    def addmeasurementcircuit(self, rounds, check='p'):
        for r in range(rounds):
            if check=='p':
                for i in self.plaquettes:
                    ancillaind = self.plaquettes[i].getancillaind()
                    self.circuit.append_operation("H", self.sqbitindmap[ancillaind])
                    self.addNoise([self.sqbitindmap[ancillaind]])
                    dqbitinds = self.plaquettes[i].getqbitinds()
                    for dind in dqbitinds:
                        self.circuit.append_operation("CNOT", [self.sqbitindmap[ancillaind], self.sqbitindmap[dind]])
                        self.addNoise([self.sqbitindmap[ancillaind], self.sqbitindmap[dind]] )
                    self.circuit.append_operation("H", self.sqbitindmap[ancillaind])
                    self.addNoise([self.sqbitindmap[ancillaind]])
            if check=='v':
                for i in self.vertices:
                    ancillaind = self.vertices[i].getancillaind()
                    dqbitinds = self.vertices[i].getqbitinds()
                    for dind in dqbitinds:
                        self.circuit.append_operation("CNOT", [ self.sqbitindmap[dind], self.sqbitindmap[ancillaind]])
                        self.addNoise([self.sqbitindmap[ancillaind], self.sqbitindmap[dind]] )
            ancillainds = sorted(self.ancillaqubits.keys())
            self.circuit.append_operation("MR", ancillainds)
            ancillalen = len(ancillainds) 
            for l in range(len(ancillainds)):
                self.circuit.append_operation("DETECTOR", [stim.target_rec(-1-l), stim.target_rec(-1-l-ancillalen)], [self.ancillaqubits[ancillainds[-1-l]].gety(), self.ancillaqubits[ancillainds[-1-l]].getx()])
        if check=='v':
            qbitinds = [self.sqbitindmap[qind] for qind in self.vertices[0].getqbitinds()]
            self.circuit.append_operation("MR", qbitinds)
        elif check=='p':
            qbitinds = [self.sqbitindmap[qind] for qind in self.plaquettes[0].getqbitinds()]
            self.circuit.append_operation("MX", qbitinds)
        self.circuit.append_operation("OBSERVABLE_INCLUDE", [stim.target_rec(-1), stim.target_rec(-2), stim.target_rec(-3), stim.target_rec(-4)], 0)
        self.circuit.append_operation("TICK")

    def addnoise(self, noise):
        if not self.initialized:
            print("ERROR topology has not been initialized")
            return
        self.circuit.append_operation("DEPOLARIZE1", self.qubits.keys(), noise)

    def enablespecialindmap(self):
        qinds = list(self.qubits.keys())
        for i in range(len(qinds)):
            for j in range(i+1,len(qinds),1):
                ivalx = self.qubits[qinds[i]].getx()
                jvalx = self.qubits[qinds[j]].getx()
                ivaly = self.qubits[qinds[i]].gety()
                jvaly = self.qubits[qinds[j]].gety()

                if ivalx > jvalx:
                    t = qinds[i]
                    qinds[i] = qinds[j]
                    qinds[j] = t
                elif ivalx == jvalx:  
                    if ivaly > jvaly:
                        t = qinds[i]
                        qinds[i] = qinds[j]
                        qinds[j] = t
        for ind,v in enumerate(qinds):
            self.sqbitindmap[v] = ind

    def parseDEM(self):
        model = self.circuit.detector_error_model(decompose_errors=True)
        return pymatching.Matching.from_detector_error_model(model)

def collect_and_plot_data(probs, error_type, rounds, shots=1000):
    if error_type == 'flip':
        checktype = 'p' 
    elif error_type == 'phase':
        checktype = 'v' 
    lerrors = {}
    numqubits = {}
    for l in [5, 8, 16, 20]:
        lerrors[l] = []
        for p in probs:
            t = Topology(l*l, l, l)
            t.setNoise(p)
            t.implementcircuit()
            t.addancilla(check=checktype)
            t.addmeasurementcircuit(rounds, check=checktype)
            numqubits[l] = t.getNumQubits()
            m = t.parseDEM()
            samples, obs = t.takeameasurement(detector=True,num_shots=shots)
            nerrors = 0
            for ind in range(shots):
                pred_ob = m.decode(samples[ind,:])
                nerrors += not np.array_equal(obs[ind, :], pred_ob)
            lerrors[l].append(nerrors/(1.0*shots))

    for l in [5, 8, 16, 20]:
        plt.loglog(probs, lerrors[l], label="qubits#={}".format(numqubits[l]))
    plt.legend()
    plt.ylabel("Logical Error Rate")
    plt.xlabel("Physical Error Rate")
    plt.title("LogLog Logical vs Physical Error Rate. Toric Code")
    plt.show()


            
if __name__ == '__main__':
    probs = [0.001, 0.01, 0.04, 0.06, 0.08, 0.1, 0.14, 0.16,  0.18, 0.2, 0.22, 0.28, 0.3 ]
    collect_and_plot_data(probs, 'flip', 2, shots=1000) 


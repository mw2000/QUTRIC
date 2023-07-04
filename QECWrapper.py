import stim
import sys
class Qubit(object):
    def __init__(self, x, y, ind):
        self.x = x
        self.y = y
        self.ind = ind
        self.processed = False

    def getind(self):
        return self.ind

    def isprocessed(self):
        return self.processed
    
    def setprocessed(self):
        self.processed = True

    def gety(self):
        return self.y

    def getx(self):
        return self.x

    def printme(self, pref=""):
        print(F'{pref} {self.x} {self.y} {self.ind} {self.processed}') 

class Plaquette(object):
    def __init__(self, locx, locy):
        self.qubitinds = []
        self.x = locx
        self.y = locy

    def addqbitinds(self, qbitinds):
        self.qubitinds = qbitinds

    def getqbitinds(self):
        return self.qubitinds

    def printme(self, pref=""):
        print(F'{pref} {self.x} {self.y}')

class Topology(object):
    def __init__(self, k, r, c):
        self.nplaquettes = k
        self.rows = r
        self.columns = c
        self.plaquettes = {}
        self.qinddict = {}
        self.qubits = {}
        self.sqbitindmap = {}
        self.qind = 0
        self.initialized = False
        self.initializetopology()
        self.circuit = stim.Circuit()

    def getqbits(self, i, j):
        xs = [i - 1, i, i + 1, i]
        ys = [j, j - 1, j, j + 1]
        rqubitinds = []
        for xbar, ybar in zip(xs, ys):
            x = xbar
            y = ybar
            if xbar % (2*self.rows) == 0:
                x = 0 
            if ybar % (2*self.columns) == 0:
                y = 0 
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
        self.initialized = True

    def implementcircuit(self):
        if not self.initialized:
            print("ERROR topology has not been initialized")
            return
        lines = []
        for i in range(self.nplaquettes-1):
            qbitinds = self.plaquettes[i].getqbitinds()
            cqbitind = -1
            for qi in qbitinds:
                if self.qubits[qi].isprocessed() == False:
                    cqbitind = qi
                    break
            qbitinds.remove(cqbitind)
            self.qubits[cqbitind].setprocessed()
            lines.append("H {}".format(self.sqbitindmap[cqbitind]))
            self.circuit.append_operation("H", self.sqbitindmap[cqbitind])
            for qi in qbitinds:
                lines.append("CNOT {} {}".format(self.sqbitindmap[cqbitind], self.sqbitindmap[qi]))
                self.circuit.append_operation("CNOT", [self.sqbitindmap[cqbitind], self.sqbitindmap[qi]])
                self.qubits[qi].setprocessed()

    def printcircuit(self):
        if not self.initialized:
            print("ERROR topology has not been initialized")
            return
        print(self.circuit)

    def addmeasurement(self):
        if not self.initialized:
            print("ERROR topology has not been initialized")
            return
        self.circuit.append_operation("MR", self.qubits.keys())
        
    def takeameasurement(self):
        if not self.initialized:
            print("ERROR topology has not been initialized")
            return
        sample = self.circuit.compile_sampler().sample(1)[0]
        print("".join(str(int(e)) for e in sample))

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
            



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    t = Topology(6, 2, 3)
    #t.addmeasurement()
    #t.takeameasurement()
    t.enablespecialindmap()
    t.implementcircuit()
    t.addmeasurement()
    #t.printcircuit()
    t.takeameasurement()


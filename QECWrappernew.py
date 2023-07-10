import stim
import sys
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

    def gety(self):
        return self.y

    def getx(self):
        return self.x

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

    def getvertexqbits(self, xv, yv):
        xs = [xv-1, xv, xv, xv+1]
        ys = [yv, yv-1, yv, yv+1]
        rqubitinds = []
        for xbar, ybar in zip(xs, ys):
            x = xbar
            y = ybar
            if xbar < 0 or xbar >= (2*self.rows) :
                x = 0 
            if ybar < 0 or ybar >= (2*self.columns):
                y = 0 
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
                x = xbar
                y = ybar
                if xbar % (2*self.rows) == 0:
                    x = 0 
                if ybar % (2*self.columns) == 0:
                    y = 0 
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

    def addancilla(self):
        self.ancillaqubits = {}
        for i in range(self.nplaquettes):
            xi = self.plaquettes[i].getx()
            yi = self.plaquettes[i].gety()
            qbit = Qubit(xi, yi, self.qind)
            qbit.setancilla("p")
            self.qinddict["{}-{}".format(xi,yi)] = self.qind 
            self.ancillaqubits[self.qind] = qbit
            self.sqbitindmap[self.qind] = self.qind
            self.plaquettes[i].addancillaqbitind(self.qind)
            self.qind += 1
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

    def addmeasurementcircuit(self):
        for i in range(self.nplaquettes):
            ancillaind = self.plaquettes[i].getancillaind()
            self.circuit.append_operation("H", self.sqbitindmap[ancillaind])
            dqbitinds = self.plaquettes[i].getqbitinds()
            for dind in dqbitinds:
                self.circuit.append_operation("CNOT", [self.sqbitindmap[ancillaind], self.sqbitindmap[dind]])
            self.circuit.append_operation("H", self.sqbitindmap[ancillaind])
        for i in self.vertices:
            ancillaind = self.vertices[i].getancillaind()
            #self.circuit.append_operation("H", self.sqbitindmap[ancillaind])
            dqbitinds = self.vertices[i].getqbitinds()
            for dind in dqbitinds:
                self.circuit.append_operation("CNOT", [ self.sqbitindmap[dind], self.sqbitindmap[ancillaind]])
        self.circuit.append_operation("MR", self.ancillaqubits.keys())

    def addnoise(self, noise):
        if not self.initialized:
            print("ERROR topology has not been initialized")
            return
        self.circuit.append_operation("X_ERROR", self.qubits.keys(), noise)

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
    #t.enablespecialindmap()
    t.implementcircuit()
    t.addancilla()
    t.addnoise(0.01)
    t.addmeasurementcircuit()
    #t.addmeasurement()
    #t.printcircuit()
    t.takeameasurement()
    print(len(t.qubits.keys()))

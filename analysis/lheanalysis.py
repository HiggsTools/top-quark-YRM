import ROOT, sys

inf = ROOT.TFile(sys.argv[1])
tree = inf.Get("events")

tree.Print("ALL")

pdgIdMap = {
    1: "d",
    2: "u",
    3: "s",
    4: "c",
    5: "b",
    6: "t",
    25: "h",
    24: "W",
    23: "Z",
    22: "gamma",
    21: "g",
    11: "e",
    12: "nue",
    13: "mu",
    14: "numu",
    15: "tau",
    16: "nutau",
}

def getParticleKind(id):
    aid = abs(id)
    r = pdgIdMap[aid]

    if r in ["W" ,"e", "mu", "tau"]:
        r = r + ("-" if id<0 else "+")
    if id < 0:
        if r in ["u", "d", "c", "s", "t", "b"]:
            r = "~" + r
    return r

class Particle:
    def __init__(self, px, py, pz, e, status, id):
        self.px = px
        self.py = py
        self.pz = pz
        self.e = e
        self.status = status
        self.id = id
        self.mothers = []

    def __str__(self):
        return "p[{6}] ({0:.2f}, {1:.2f}, {2:.2f}, {3:.2f}) {4} {5}".format(
            self.px, self.py, self.pz, self.e, self.status,
            [getParticleKind(m.id) for m in self.mothers ],
            getParticleKind(self.id)
        )

nb = 0
nproc = 0
for iev in range(tree.GetEntries()):
    nb += tree.GetEntry(iev)
    nparticles = tree.n_particles

    particles = []
    particles_s1 = []
    particles_s2 = []

    for ipt in range(nparticles):
        pid = tree.PID[ipt]
        px = tree.P_X[ipt]
        py = tree.P_Y[ipt]
        pz = tree.P_Z[ipt]
        e = tree.E[ipt]
        mid1 = tree.MID1[ipt]
        mid2 = tree.MID2[ipt]
        status = tree.STATUS[ipt]

        part = Particle(px, py, pz, e, status, pid)
        particles += [part]
        part.mothers = [mid1]
        if mid2!=mid1:
            part.mothers += [mid2]
    
    for part in particles:
        mids = part.mothers
        #Assign mothers
        #if mid = 0 -> no mother in LHE file
        #lhe indexing 1-based
        part.mothers = [particles[mid-1] for mid in mids if mid-1>=0]

        if part.status == 1:
            particles_s1 += [part]
        elif part.status == 2:
            particles_s2 += [part]
    
    print "---"
    for p in particles:
        print str(p)
    nproc += 1

print "processed {0} entries, read {1:.2f} Mb".format(nproc, nb/1024.0/1024.0)

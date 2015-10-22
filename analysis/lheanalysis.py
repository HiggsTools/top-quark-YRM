import ROOT, sys
import math

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

out=ROOT.TFile("analysis.root","recreate")


#Histograms
Cosl1l2=ROOT.TH1D("Cosl1l2", "CosTheta", 20,-1,1)



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
    
    #print "---"
    l1=ROOT.TLorentzVector()
    l2=ROOT.TLorentzVector()
    gamma1=ROOT.TLorentzVector()
    gamma2=ROOT.TLorentzVector()


    for p in particles:
        if abs(p.id) is 22:
            if gamma1.Mag() !=0:
                gamma1.SetPx(p.px)
                gamma1.SetPy(p.py)
                gamma1.SetPz(p.pz)
                gamma1.SetE(p.e)
                if gamma1.Pt() <= 20:
                    continue
                if gamma1.Eta() > 2.5:
                    continue

        else:
                gamma2.SetPx(p.px)
                gamma2.SetPy(p.py)
                gamma2.SetPz(p.pz)
                gamma2.SetE(p.e)
                if gamma2.Pt() <= 20:
                    continue
                if gamma2.Eta() > 2.5:
                    continue
            
        
        
        if abs(p.id) in [11,13]:
            mothers=p.mothers
            for mother in mothers:
                if 24 is abs(mother.id):
                    Wmothers=mother.mothers
                    for Wmother in Wmothers:
                        if 6 is abs(Wmother.id):
                            l=ROOT.TLorentzVector()
                            l.SetPx(p.px)
                            l.SetPy(p.py)
                            l.SetPz(p.pz)
                            l.SetE(p.e)
                            
                            Mother=ROOT.TLorentzVector()
                            Mother.SetPx(Wmother.px)
                            Mother.SetPy(Wmother.py)
                            Mother.SetPz(Wmother.pz)
                            Mother.SetE(Wmother.e)

                            l.Boost(Mother.BoostVector())
                            if p.id > 0:
                                l1=l
                            else:
                                l2=l


    

    if gamma1.DeltaR(gamma2) < 0.4:
        continue
    if (gamma1+gamma2).M() < 123 and (gamma1+gamma2).M() > 129:
        continue
    if l1.DeltaR(l2) < 0.4:
        continue
    #if l1.Pt() < 5 or l2.Pt() < 5:
    #    continue



    Cosl1l2.Fill(math.cos(l1.Angle(l2.Vect())))


    nproc += 1
out.Write()

print "processed {0} entries, read {1:.2f} Mb".format(nproc, nb/1024.0/1024.0)

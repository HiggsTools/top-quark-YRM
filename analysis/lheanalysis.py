import ROOT, sys
import math
import numpy as np
import copy

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

out=ROOT.TFile(sys.argv[2],"recreate")


#Histograms
Cosl1l2=ROOT.TH1D("Cosl1l2", "CosTheta", 20,-1,1)
Cosl1l2_frame2=ROOT.TH1D("Cosl1l2Frame2", "CosTheta", 20,-1,1)

Etal_LabFrame=ROOT.TH1D("Etal_LabFrame","#Eta_{l} in LabFrame", 24, 0,3)
Etab_LabFrame=ROOT.TH1D("Etab_LabFrame","#Eta_{b} in LabFrame", 24, 0,3)
Cosl1l2_LabFrame=ROOT.TH1D("Cosl1l2_LabFrame","#Cos_{ll} in LabFrame", 24, 0,3)
Cosb1b2_LabFrame=ROOT.TH1D("Cosb1b2_LabFrame","#Cos_{bb} in LabFrame", 24, 0,3)


hleps = ROOT.TH1D("LepId", "LepId", 20, 0, 20)
hlep1pt = ROOT.TH1D("Lep1Pt", "First lepton pt", 20, 0, 300)
hlep2pt = ROOT.TH1D("Lep2Pt", "second lepton pt", 20, 0, 300)
hlepdr = ROOT.TH1D("LepDr", "lepton dr", 20, 0, 5)


outree = ROOT.TTree("events", "events")
arr_cosl1l2 = np.zeros(1, "d")
arr_cosl1l2_frame2 = np.zeros(1, "d")
arr_lep1_pt = np.zeros(1, "d")
arr_lep2_pt = np.zeros(1, "d")
arr_top1_pt = np.zeros(1, "d")
arr_top2_pt = np.zeros(1, "d")
arr_lep1_pt_frame2 = np.zeros(1, "d")
arr_lep2_pt_frame2 = np.zeros(1, "d")

outree.Branch("cosl1l2", arr_cosl1l2, "cosl1l2/D")
outree.Branch("cosl1l2_frame2", arr_cosl1l2_frame2, "cosl1l2_frame2/D")
outree.Branch("lep1_pt_frame2", arr_lep1_pt_frame2, "lep1_pt_frame2/D")
outree.Branch("lep2_pt_frame2", arr_lep2_pt_frame2, "lep2_pt_frame2/D")
outree.Branch("lep1_pt", arr_lep1_pt, "lep1_pt/D")
outree.Branch("lep2_pt", arr_lep2_pt, "lep2_pt/D")
outree.Branch("top1_pt", arr_top1_pt, "top1_pt/D")
outree.Branch("top2_pt", arr_top2_pt, "top2_pt/D")

#for iev in range(tree.GetEntries()):
for iev in range(50000):
    if iev%1000==0:
        print "event", iev
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
    b1=ROOT.TLorentzVector()
    b2=ROOT.TLorentzVector()

    

    top1 = None
    top2 = None

    for np, p in enumerate(particles):
        if abs(p.id) is 22:
            if gamma1.Mag() == 0:
                gamma1.SetPx(p.px)
                gamma1.SetPy(p.py)
                gamma1.SetPz(p.pz)
                gamma1.SetE(p.e)
                if gamma1.Pt() <= 20 or gamma1.Eta() > 2.5:
                    break
            else:
                gamma2.SetPx(p.px)
                gamma2.SetPy(p.py)
                gamma2.SetPz(p.pz)
                gamma2.SetE(p.e)
                if gamma2.Pt() <= 20 or gamma2.Eta() > 2.5:
                    break
        
        
        if abs(p.id) in [11, 13]:
            hleps.Fill(p.id)
            #print np, p.id
            mothers=p.mothers
            for mother in mothers:
                if 24 is abs(mother.id):
                    #print "found W as mother"
                    Wmothers=mother.mothers
                    for Wmother in Wmothers:
                        if 6 is abs(Wmother.id):
                            #print "found top as mother"``
                            #lepton
                            l=ROOT.TLorentzVector()
                            l.SetPx(p.px)
                            l.SetPy(p.py)
                            l.SetPz(p.pz)
                            l.SetE(p.e)

                            #top that is the mother of lepton
                            Mother=ROOT.TLorentzVector()
                            Mother.SetPx(Wmother.px)
                            Mother.SetPy(Wmother.py)
                            Mother.SetPz(Wmother.pz)
                            Mother.SetE(Wmother.e)

                            if p.id > 0:
                                l1=copy.deepcopy(l)
                                top1 = Mother
                                arr_lep1_pt[0] = l1.Pt()
                                arr_top1_pt[0] = Mother.Pt()
                            else:
                                l2=copy.deepcopy(l)
                                top2 = Mother
                                arr_lep2_pt[0] = l2.Pt()
                                arr_top2_pt[0] = Mother.Pt()
        if 5 is abs(p.id):
            mothers=p.mothers
            for mother in mothers:
                    if 6 is abs(mother.id):
                        b=ROOT.TLorentzVector()
                        b.SetPx(p.px)
                        b.SetPy(p.py)
                        b.SetPz(p.pz)
                        b.SetE(p.e)
                   
                        if p.id > 0:
                            b1=copy.deepcopy(b)
                        else:
                            b2=copy.deepcopy(b)


    if gamma1.Mag() == 0 or gamma2.Mag() == 0 or l1.Mag()==0 or l2.Mag()==0 or b1.Mag() == 0 or b2.Mag() == 0:
        continue #next event

    toppair = top1 + top2
    
    #print (l1-l2).Mag()

    l1_frame1 = ROOT.TLorentzVector(l1)
    l2_frame1 = ROOT.TLorentzVector(l2)

    l1_frame1.Boost(-toppair.BoostVector())
    l2_frame1.Boost(-toppair.BoostVector())

    l1_frame1.Boost(-top1.BoostVector())
    l2_frame1.Boost(-top2.BoostVector())

    l1_frame2 = ROOT.TLorentzVector(l1)
    l1_frame2.Boost(-top1.BoostVector())

    l2_frame2 = ROOT.TLorentzVector(l2)
    l2_frame2.Boost(-top2.BoostVector())
    
    b1_LabFrame = ROOT.TLorentzVector(b1)
    b2_LabFrame = ROOT.TLorentzVector(b2)

    l1_LabFrame = ROOT.TLorentzVector(l1)
    l2_LabFrame = ROOT.TLorentzVector(l2)

    # if gamma1.DeltaR(gamma2) < 0.4:
    #     continue
    # if (gamma1+gamma2).M() < 123 or (gamma1+gamma2).M() > 129:
    #     continue
    # if l1.DeltaR(l2) < 0.4:
    #     continue
    if l1.Pt() < 10 or l2.Pt() < 10:
       continue
    if abs(l1.Eta()) > 2.5 or abs(l2.Eta())>2.5:
        continue
    hlep1pt.Fill(l1.Pt())
    hlep2pt.Fill(l2.Pt())

    hlepdr.Fill(l1.DeltaR(l2))

    #print "l1", l1.Pt(), l1.Eta(), l1.Phi(), l1.M()
    #print "l2", l2.Pt(), l2.Eta(), l2.Phi(), l2.M()

    cosl1l2 = math.cos(l1_frame1.Angle(l2_frame1.Vect()))
    Cosl1l2.Fill(cosl1l2)

    cosl1l2_frame2 = math.cos(l1_frame2.Angle(l2_frame2.Vect()))
    Cosl1l2_frame2.Fill(cosl1l2_frame2)
    
    Etal_LabFrame.Fill(abs(l1.Eta()-l2.Eta()))
    Etab_LabFrame.Fill(abs(b1.Eta()-b2.Eta()))
   
    cosl1l2_LabFrame = math.cos(l1_LabFrame.Angle(l2_LabFrame.Vect()))
    Cosl1l2_LabFrame.Fill(cosl1l2_LabFrame)

    cosb1b2_LabFrame = math.cos(b1_LabFrame.Angle(b2_LabFrame.Vect()))
    Cosb1b2_LabFrame.Fill(cosb1b2_LabFrame)


    #print cosl1l2, cosl1l2_frame2

    arr_cosl1l2[0] = cosl1l2
    arr_cosl1l2_frame2[0] = cosl1l2_frame2

    arr_lep1_pt_frame2[0] = l1.Pt()
    arr_lep2_pt_frame2[0] = l2.Pt()
    outree.Fill()

    nproc += 1
out.Write()


print "processed {0} entries, read {1:.2f} Mb".format(nproc, nb/1024.0/1024.0)

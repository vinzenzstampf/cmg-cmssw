from PhysicsTools.Heppy.physicsobjects.PhysicsObjects import printOut 
from PhysicsTools.Heppy.physicsobjects.PhysicsObjects import GenParticle 
import ROOT

def findStatus1Leptons(particle):
    '''Returns status 1 e and mu among the particle daughters'''
    leptons = []
    for i in range( particle.numberOfDaughters() ):
        dau = particle.daughter(i)
        if dau.status() == 1:
            if abs(dau.pdgId())==11 or abs(dau.pdgId())==13:
                leptons.append( dau )
            else:
                continue
        else:
            leptons = findStatus1Leptons( dau, leptons )
    return leptons


def allDaughters(particle, daughters, rank ):
    '''Fills daughters with all the daughters of particle.
    Recursive function.'''
    rank += 1 
    for i in range( particle.numberOfDaughters() ):
        dau = GenParticle(particle.daughter(i))
        dau.rank = rank
        daughters.append( dau )
        daughters = allDaughters( dau, daughters, rank )
    return daughters


def bosonToX(particles, bosonType, xType):
    bosons = [x for x in particles if x.status()==3 and x.pdgId()==bosonType]
    daughters = []
    if len(bosons)==0:
        return [], False
    boson = bosons[0]
    daus = []
    allDaughters( boson, daus, 0)
    xDaus = [x for x in daus if x.status()==3 and abs(x.pdgId())==xType]
    # print printOut(xDaus)
    return xDaus, True 

def isNotHadronicId(pdgId,includeSMLeptons=True):
    if abs(pdgId) in [11,12,13,14,15,16]:
        return includeSMLeptons
    i = (abs(pdgId) % 1000)
    return i > 10 and i != 21 and i < 100

def isPromptLepton(lepton, beforeFSR, includeMotherless=True, includeTauDecays=False):
    if abs(lepton.pdgId()) not in [11,13,15]:
        return False
    if lepton.numberOfMothers() == 0:
        return includeMotherless;
    mom = lepton.mother()
    if mom.pdgId() == lepton.pdgId():
        if beforeFSR: return False
        return isPromptLepton(mom, beforeFSR, includeMotherless, includeTauDecays)
    elif abs(mom.pdgId()) == 15:
        if not includeTauDecays: return False
        return isPromptLepton(mom, beforeFSR, includeMotherless, includeTauDecays)
    else:
        return isNotHadronicId(mom.pdgId(), includeSMLeptons=False)


def isNotFromHadronicShower(l):
    for x in xrange(l.numberOfMothers()):
        mom = l.mother(x)
        if mom.status() > 2: return True
        id = abs(mom.pdgId())
        if id > 1000000: return True
        if id > 100: return False
        if id <   6: return False
        if id == 21: return False
        if id in [11,12,13,14,15,16]: 
            if l.status() > 2: return True
            return isNotFromHadronicShower(mom)
        if id >= 22 and id <= 39: return True
    return True

def realGenDaughters(gp,dauChecked=[],excludeRadiation=True):
    """Get the daughters of a particle, going through radiative X -> X' + a
       decays, either including or excluding the radiation among the daughters
       e.g. for  
                  X -> X' + a, X' -> b c 
           realGenDaughters(X, excludeRadiation=True)  = { b, c }
           realGenDaughters(X, excludeRadiation=False) = { a, b, c }
       Keep track of the daughters we've already checked so that we don't get
       stuck in an infinite loop"""
    ret = []
    for i in xrange(gp.numberOfDaughters()):
        dau = gp.daughter(i)
        if dau in dauChecked: continue
        dauChecked.append(dau)
        if dau.pdgId() == gp.pdgId():
            if excludeRadiation:
                return realGenDaughters(dau, dauChecked)
            else:
                ret += realGenDaughters(dau, dauChecked)
        else:
            ret.append(dau)
    return ret

def realGenMothers(gp):
    """Get the mothers of a particle X going through intermediate X -> X' chains.
       e.g. if Y -> X, X -> X' realGenMothers(X') = Y"""
    ret = []
    for i in xrange(gp.numberOfMothers()):
        mom = gp.mother(i)
        if mom.pdgId() == gp.pdgId():
            ret += realGenMothers(mom)
        else:
            ret.append(mom)
    return ret

def motherRef(gp,i=0):
    return (gp.mother(), ROOT.heppy.GenParticleRefHelper.motherKey(gp,i))
def daughterRef(gp,i=0):
    return (gp.daughter(), ROOT.heppy.GenParticleRefHelper.daughterKey(gp,i))

def lastGenCopy(gp):
    me = gp.pdgId();
    for i in xrange(gp.numberOfDaughters()):
        if gp.daughter(i).pdgId() == me:
            return False
    return True



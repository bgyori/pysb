from pysb import *
import pysb.generator.bng as bng

# from http://bionetgen.org/index.php/Simple_EGFR_model


Model('simple_egfr')



# Concentrations in number per cell
Parameter('EGF_tot',      1.2e6)
Parameter('EGFR_tot',     1.8e5)
Parameter('Grb2_tot',     1.5e5)
Parameter('Sos_tot',      6.2e4)

# Biomolecular rate constants are in (# per cell)^-1 s^-1,
#  obtained by dividing constants in M^-1 s^-1 by Na*V,
#  where Na is Avogadro's number and V is the volume
#  of the relevant compartment (the cytoplasm for all cases here).
# Unimolecular rate constants are in s^-1
Parameter('kp1',      1.667e-06) # ligand-monomer binding
Parameter('km1',           0.06) # ligand-monomer dissociation

Parameter('kp2',      5.556e-06) # aggregation of bound monomers
Parameter('km2',            0.1) # dissociation of bound monomers

Parameter('kp3',            0.5) # dimer transphosphorylation   
Parameter('km3',          4.505) # dimer dephosphorylation        

Parameter('kp4',      8.333e-07) # binding of Grb2 to receptor
Parameter('km4',           0.05) # dissociation of Grb2 from receptor

Parameter('kp5',      5.556e-06) # binding of Grb2 to Sos
Parameter('km5',           0.06) # dissociation of Grb2 from Sos

Parameter('kdeg',          0.01)



#  EGF(r)
Monomer('EGF', 'r')

#  EGFR(l, r, Y1068~U~P, Y1148~U~P)
Monomer('EGFR',
        ['l','r','Y1068','Y1148'],
        { 'Y1068': ['U','P'],
          'Y1148': ['U','P'] }
        )

#  Grb2(SH2, SH3)
Monomer('Grb2', ['SH2','SH3'])

#  Shc(PTB, Y317~U~P)
Monomer('Shc', ['PTB','Y317'], { 'Y317': ['U','P'] } )

#  Sos(PR)
Monomer('Sos', 'PR')

#  NULL()
Monomer('NULL')



# Ligand-receptor binding      
#EGFR(l,r) + EGF(r) <-> EGFR(l!1,r).EGF(r!1) kp1, km1
Rule('ligand_receptor_binding',
     EGFR(l=None, r=None) + EGF(r=None) <>
     EGFR(l=1, r=None)    * EGF(r=1),
     kp1, km1)

# Receptor-aggregation 
#EGFR(l!+,r) + EGFR(l!+,r) <-> EGFR(l!+,r!1).EGFR(l!+,r!1) kp2,km2
Rule('receptor_aggregation',
     EGFR(l=ANY, r=None) + EGFR(l=ANY, r=None) <>
     EGFR(l=ANY, r=1)    * EGFR(l=ANY, r=1),
     kp2, km2)

# Transphosphorylation of EGFR by RTK
#EGFR(r!+,Y1068~U) -> EGFR(r!+,Y1068~P)  kp3
Rule('transphos_egfr',
     EGFR(r=ANY, Y1068='U') >>
     EGFR(r=ANY, Y1068='P'),
     kp3)

# Dephosphorylayion
#EGFR(Y1068~P) -> EGFR(Y1068~U)  km3
Rule('dephos_egfr',
     EGFR(Y1068='P') >>
     EGFR(Y1068='U'),
     km3)

# Grb2 binding to pY1068
#EGFR(Y1068~P) + Grb2(SH2)   <-> EGFR(Y1068~P!1).Grb2(SH2!1)   kp4,km4
Rule('grb2_bind_egfr',
     EGFR(Y1068='P')     + Grb2(SH2=None) <>
     EGFR(Y1068=('P',1)) * Grb2(SH2=1),
     kp4, km4)

# Grb2 binding to Sos
#Grb2(SH2,SH3) + Sos(PR) <-> Grb2(SH2,SH3!1).Sos(PR!1) kp5,km5
Rule('grb2_bind_sos',
     Grb2(SH2=None, SH3=None) + Sos(PR=None) <>
     Grb2(SH2=None, SH3=1)    * Sos(PR=1),
     kp5, km5)


#  D        EGFR(l!+)
simple_egfr.observe('D', EGFR(l=ANY))
#  RP	   EGFR(Y1068~P!?)
simple_egfr.observe('RP', EGFR(Y1068=('P',WILD)))
#  R_Grb2   EGFR(Y1068!1).Grb2(SH2!1)
simple_egfr.observe('R_Grb2', EGFR(Y1068=1) * Grb2(SH2=1))
#  Sos_act  EGFR(Y1068!1).Grb2(SH2!1,SH3!2).Sos(PR!2)
simple_egfr.observe('Sos_act', EGFR(Y1068=1) * Grb2(SH2=1, SH3=2) * Sos(PR=2))
#  EGFR_tot EGFR()
simple_egfr.observe('EGFR_tot', EGFR())


gen = bng.BngGenerator(model=simple_egfr)
print gen.content
print """
begin species
EGF(r)                           EGF_tot
EGFR(l,r,Y1068~U,Y1148~U)        EGFR_tot
Grb2(SH2,SH3)                    Grb2_tot
Sos(PR)                          Sos_tot
end species

begin actions
generate_network({overwrite=>1});
end actions
"""

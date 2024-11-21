# Julia code for interaction-mediated entanglement in
# Continuous-Variable gaussian states
# Copyright (C) 2024  Ankit Kumar
# Email: kumar.ankit.vyas@gmail.com
#-----------------------------------------------------------------------
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#-----------------------------------------------------------------------


#Global Parameters
#-------------------------------------------------- 
setprecision(256)
#-------------------------------------------------- 

const hb = 1.05457182e-34
const c = 2.99792458e8
const G = 6.6743e-11
const kB = 1.380649e-2

const rho_Os = 22.59e3
const rho_Pt = 21.45e3
const rho_SiO2 = 2.65e3


# calculate mass from radius
#--------------------------------------------------
function Mass_f_Rad(rad,rho)
    return rho * 4/3*pi*rad^3
end


# calculate radius from mass
#--------------------------------------------------
function Rad_f_Mass(mass,rho)
    return ( 3*mass/(4*pi*rho) )^(1.0/3)
end


# calculate entanglement at a given time
# outputs = Logarithmic Negativity, Entanglement Entropy
#--------------------------------------------------
function Calc_Entanglement(t,m,w0,w,nbar=0)

    detA,detB,detG,detCM = CovMat_FreeFall(t,m,w0,w)
    LogNeg,Entropy = Entanglement_CovMat(detA,detB,detG,detCM,nbar)

    #-------------------------------
    return LogNeg,Entropy
end


# retrn values of entanglement given a Covariance matrix
#--------------------------------------------------
function Entanglement_CovMat(detA,detB,detG,detCM,nbar=0)
    
    SumT = detA+detB-2*detG    
    numT = sqrt( SumT - sqrt(SumT^2-4*detCM) ) / sqrt(2)
    LogNeg = max( 0, -log2( (2*nbar+1)*2*numT ) )
    
    ee = (2*nbar+1)*sqrt(detA)
    Entropy = (ee+1/2)*log2(ee+1/2) - (ee-1/2)*log2(ee-1/2)

    #-------------------------------
    return LogNeg,Entropy
end


# returns the covariance matrix of two interacting masses in free fall
# expressions defined in https://doi.org/10.22331/q-2023-05-15-1008
#--------------------------------------------------
function CovMat_FreeFall(t,m,w0,w)
    
    a = 1.0/(4*m*w0) * ( 2 + (w0*t)^2 + (1+(w0/w)^2 ) * sinh(w*t)^2 )   # sig_00
    b = 1.0/(4*m*w0) * ( (w0*t)^2 - (1+(w0/w)^2 ) * sinh(w*t)^2 )       # sig_02
    c = (m*w0/4) * ( 2 + (1+(w/w0)^2 ) * sinh(w*t)^2 )                  # sig_11
    d = -(m*w0/4) * (1+(w/w0)^2) * sinh(w*t)^2                          # sig_13
    e = 1.0/8 *( 2*w0*t + (w0/w+w/w0) * sinh(2*w*t) )                   # sig_01
    f = 1.0/8 *( 2*w0*t - (w0/w+w/w0) * sinh(2*w*t) )                   # sig_03

    detA = a*c-e^2          # | alpha |
    detB = detA             # | beta |
    detG = b*d-f^2          # | gamma |
    detCM = 1.0/16          # | Covariance Matrix sigma |

    #-------------------------------
    return detA,detB,detG,detCM
end


# entanglement quantifier for the Casimir interaction
# wC_sq is defined in https://doi.org/10.1103/PhysRevD.109.L101501
# form of interaction and constants taken from https://doi.org/10.1103/PhysRevLett.99.170403
#--------------------------------------------------
function get_wC_sq(R,L,m)

    n_terms = 10
    cc = zeros(n_terms)
    cc[1] = 143. / 16
    cc[2] = 0
    cc[3] = 7947. / 160
    cc[4] = 2065. / 32
    cc[5] = 27705347. / 100800.
    cc[6] = -55251. / 64.
    cc[7] = 1373212550401. / 144506880.
    cc[8] = -7583389. / 320.
    cc[9] = -2516749144274023. / 44508119040.
    cc[10] = 274953589659739. / 275251200.

    wC_sq = 0
    for n in range(0,n_terms-1)
        wC_sq += cc[n+1]*(n^2+15*n+56)*(R/L)^n  
    end
    wC_sq *= 2*hb*c*R^6/(m*pi*L^9)

    #-------------------------------
    return wC_sq
end



#--------------------------------------------------
#CODE INPUTS
#--------------------------------------------------
rho = rho_Pt

m = 0.25*1e-12*1e-3
R0 = Rad_f_Mass(m,rho)

L = 2.5*R0
sig = 2.5e-9
w0 = hb/(2*m*sig^2)

tmax = 5
npts = 1+10

t_list = range(0,tmax,npts)


# omega^2 for different interactions
#----------------------------------

wN_sq = 4*G*m/L^3
#wM_sq = 8/3*(sqrt(2)-1)*sqrt(G*m*a0)/L^2
#wM_sq = get_wC_sq(R0,L,m)

#wNC_sq = wN_sq + wM_sq
#wMC_sq = wM_sq + wM_sq


# choose the interaction
#--------------------------------------------
w_sq = wN_sq


# calculations for given Parameters
#--------------------------------------------
w = sqrt(wN_sq)
lneg  = zeros(npts)
vnent = zeros(npts)
table = String[]
for idx in range(1,npts)
    lneg[idx],vnent[idx] = Calc_Entanglement(t_list[idx],m,w0,w)
    line = "$(t_list[idx]) \t $(lneg[idx]) \t $(vnent[idx]) \n"
    push!(table,line)
    print( line )
end


#writing the data in output file
#-------------------------------------------------------------
open("output-Entanglement.dat","w") do out_file
    write(out_file, "Time (s) \t Logarithmic Negativity \t von-Neumann Entropy \n"  )
    for line in table
        write(out_file,line)
    end

end

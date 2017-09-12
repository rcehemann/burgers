## generates a structure along the burger path between BCC and HCP
## NOTE: BCC := l1 = l2 = 0; HCP := l1 = l2 = 1.
##################################################################

import atomout
import numpy as np
from numpy import cos, sin, arccos, arctan2
from StringIO import StringIO
import sys

## function alpha returns common factor alpha(l1)
def alpha ( l1, chi ):
	return ( 1 + ( chi - 1 ) * l1 )

## function tau returns unit cell (positions divided by a_bcc) given lambda_1 and lambda_2
def tau ( l1, l2, chi, eta ):

	taucell = np.array([[0., 0., 0.],
			[0., 0., 0.],
			[0., 0., 0.],
			[0., 0., 0.]])

	taucell[0] += np.array([ 0.  , ((3. + l2)/12.) * alpha(l1, chi) * np.sqrt(2)       , (l1 * eta + np.sqrt(2))/4.      ])
	taucell[1] += np.array([ 0.  , (1. - (3. + l2)/12.) * alpha(l1, chi) * np.sqrt(2)  , 3. * (l1 * eta + np.sqrt(2))/4. ])
	taucell[2] += np.array([ 1./(2 * alpha(l1, chi)), (1./2 - (3. + l2)/12.) * alpha(l1, chi) * np.sqrt(2), 3. * (l1 * eta + np.sqrt(2))/4. ])
	taucell[3] += np.array([ 1./(2 * alpha(l1, chi)), (1./2 + (3. + l2)/12.) * alpha(l1, chi) * np.sqrt(2), (l1 * eta + np.sqrt(2))/4.      ])

	return taucell

## function taubox returns box vectors of tau unit cell
def taubox ( l1, chi, eta ):

	taubox = np.array([[0., 0., 0.],
			[0., 0., 0.],
			[0., 0., 0.]])

	taubox[0][0] += 1. / alpha(l1, chi)
	taubox[1][1] += np.sqrt(2) * alpha(l1, chi)
	taubox[2][2] += (l1 * eta + np.sqrt(2))

	return taubox

##############################################################
##		BEGIN MAIN CODE				
##############################################################

if len(sys.argv) != 7:
	print "INPUTS:"
	print "	\t 1) BCC lattice constant a"
	print "	\t 1) HCP lattice constant a"
	print "	\t 1) HCP lattice ratio  c/a"
	print "	\t 2) lambda_1 parameter"
	print "	\t 3) lambda_2 parameter"
	print " \t 4) output format, e.g. xyz"
	sys.exit()

a_bcc  = float(sys.argv[1])
a_hcp  = float(sys.argv[2])
ca_hcp = float(sys.argv[3])
lam_1  = float(sys.argv[4])
lam_2  = float(sys.argv[5])
FRM    = str  (sys.argv[6])

#chi_0  = a_bcc / a_hcp			# find notes on desk for derivation
#eta_0  = ca_hcp / chi_0 - np.sqrt(2)	#

chi_0   = np.power(3./2, 1./4)
eta_0   = 0

if lam_1 <= 1 and lam_2 <= 1:

	atoms = tau( lam_1, lam_2, chi_0, eta_0 )
	box   = taubox ( lam_1, chi_0, eta_0 )
	
	if FRM == "xyz":
		line = "SHEAR L1= " + str(lam_1) + " SHUFFLE L2= " + str(lam_2)
		atomout.output_xyz(a_bcc*atoms, line)
	elif FRM == "car":
		line = "SHEAR L1= " + str(lam_1) + " SHUFFLE L2= " + str(lam_2)
		atomout.output_car(box, a_bcc, atoms, line)
	elif FRM == "carsd":
		line = "SHEAR L1= " + str(lam_1) + " SHUFFLE L2= " + str(lam_2)
		atomout.output_car_sd(box, atoms, line, flags)
	elif FRM == "xyzcol":
		atomout.output_xyz_color(a_bcc*atoms, a_bcc*box[2][2]/2, a_bcc*cz, a_bcc*box[2][2])
	elif FRM == "lmp":
		nats = np.size(atoms, axis=0)
		line = "SHEAR L1= " + str(lam_1) + " SHUFFLE L2= " + str(lam_2)
		atomout.output_lmp(a_bcc*box, a_bcc*atoms, nats, line)
	elif FRM == "pdb":
		nats = np.size(atoms, axis=0)
		atomout.output_pdb(a_bcc*box, a_bcc*atoms, nats)
	elif FRM == "conf":
		nats = np.size(atoms, axis=0)
		atomout.output_conf(a_bcc*box, a_bcc*atoms, nats)

elif lam_1 > 1 and lam_2 > 1:
	if lam_1 == int(lam_1) and lam_2 == int(lam_2):

		for i in range(0, int(lam_1)+1):
			for j in range(0, int(lam_2)+1):

				lami = float(i)/lam_1
				lamj = float(j)/lam_2

				atoms = tau( lami, lamj, chi_0, eta_0 )
				box   = taubox( lami, chi_0, eta_0 )


				if FRM == "xyz":
					line = "SHEAR L1= " + str(lam_1) + " SHUFFLE L2= " + str(lam_2)
					atomout.output_xyz(a_bcc*atoms, line)
				elif FRM == "car":
					line = "SHEAR L1= " + str(lam_1) + " SHUFFLE L2= " + str(lam_2)
					atomout.output_car(box, a_bcc, atoms, line)
				elif FRM == "carsd":
					line = "SHEAR L1= " + str(lam_1) + " SHUFFLE L2= " + str(lam_2)
					atomout.output_car_sd(box, atoms, line, flags)
				elif FRM == "xyzcol":
					atomout.output_xyz_color(a_bcc*atoms, a_bcc*box[2][2]/2, a_bcc*cz, a_bcc*box[2][2])
				elif FRM == "lmp":
					nats = np.size(atoms, axis=0)
					line = "SHEAR L1= " + str(lam_1) + " SHUFFLE L2= " + str(lam_2)
					atomout.output_lmp(a_bcc*box, a_bcc*atoms, nats, line)
				elif FRM == "pdb":
					nats = np.size(atoms, axis=0)
					atomout.output_pdb(a_bcc*box, a_bcc*atoms, nats)
				elif FRM == "conf":
					nats = np.size(atoms, axis=0)
					atomout.output_conf(a_bcc*box, a_bcc*atoms, nats)
	else:
		print "Requires integer mesh size in each direction (args 4 and 5)"
		exit()
else:
	print "Requires lambda_1 and lambda_2 to both be less than one or integers greater than one..."
	exit()

import sys, numpy as np

# Units
# [length] = [time] = [energy]^-1 = [mass]^[-1]

# Constants
hbar 			= 1.0545718 * 10**(-34) 	# m^2 kg / s
c 				= 299792458 				# m / s
hbarc			= hbar*c 					# m^3 kg / s^2 = joules m

# hbarc = 0.19732697 #eV micro m

#### Description of electron volt:
# By definition, it is the amount of energy gained (or lost) by the charge of a
# single electron moving across an electric potential difference of one volt.
electron_volt	= 1.60217662 * 10**(-19)	# joules

# Length dictionary
length_dictionary 	= {"fermi" : 10**(-15), "nano" : 10**(-9), "micro" : 10**(-6)}

# Energy dictionary
energy_dictionary 	= {"mev" : 10**6, "gev" : 10**9}

# LENGTHS: hbarc in different unit lengths
hbarc_meters_ev		= hbarc / electron_volt 						# meters
hbarc_micro_ev		= hbarc_meters_ev * length_dictionary["micro"]	# micro-meters
hbarc_nano_ev		= hbarc_meters_ev * length_dictionary["nano"]	# nano-meters
hbarc_fermi_ev		= hbarc_meters_ev * length_dictionary["fermi"]	# fermi-meters

# ENERGY: hbarc in different energy units
hbarc_meters_mev 	= hbarc_meters_ev / energy_dictionary["mev"]
hbarc_meters_gev 	= hbarc_meters_ev / energy_dictionary["gev"]

hbarc_micro_mev		= hbarc_meters_mev / length_dictionary["micro"]
hbarc_micro_gev		= hbarc_meters_gev / length_dictionary["micro"]

hbarc_nano_mev		= hbarc_meters_mev / length_dictionary["nano"]
hbarc_nano_gev		= hbarc_meters_gev / length_dictionary["nano"]

hbarc_fermi_mev		= hbarc_meters_mev / length_dictionary["fermi"]
hbarc_fermi_gev		= hbarc_meters_gev / length_dictionary["fermi"]	# fermi-meters

one_fermi_from_mev = 1 / hbarc_fermi_mev # 1 fermi
one_fermi_from_gev = 1 / hbarc_fermi_gev # 1 fermi

one_mev_from_fermi = hbarc_fermi_mev
one_gev_from_fermi = hbarc_fermi_gev

# S = hbarc[fermi]

def getLatticeSpacing(beta):
	# if beta < 5.7: raise ValueError("Beta should be larger than 5.7!")
	r = 0.5
	bval = (beta - 6)
	a = np.exp(-1.6805 - 1.7139*bval + 0.8155*bval**2 - 0.6667*bval**3)*0.5
	return a # fermi

def convert_meters_to_ev():
	return hbarc*fm

def convert_ev_to_meters(fm):
	return fm

print "Lattice spacing: %.4f [fermi]" % getLatticeSpacing(6.45)
print "Energy scale:	%.4f [GeV]" % (getLatticeSpacing(6.45) / hbarc_fermi_gev)


def main(args):
	None
	# if args[0].lower() == "fm":
	# 	print "%g [fm] = %g [GeV]" % (float(args[1]), convert_fm_to_gev(float(args[1])))
	# elif args[0].lower() == "gev":
	# 	print "%g [GeV] = %g [fm]" % (float(args[1]), convert_gev_to_fm(float(args[1])))
	# else:
	# 	print "Arguments not recognized: ", args

if __name__ == '__main__':
	if len(sys.argv) == 1:
		args = ["fm","0.5"]
	else:
		args = sys.argv[1:]

	main(args)
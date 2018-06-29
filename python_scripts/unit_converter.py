import sys, numpy as np

# UNITS IN NATURAL UNITS

# Units
# [length] = [time] = [energy]^{-1} = [mass]^{-1}

# Constants
hbar 			= 1.0545718 * 10**(-34) 	# m^2 kg / s
c 				= 299792458 				# m / s
hbarc			= hbar*c 					# m^3 kg / s^2 = joules m

# hbarc = 0.19732697 #eV micro m

# 1 eV^{-1} =  1.97*10^{-7} m

#### Description of electron volt:
# By definition, it is the amount of energy gained (or lost) by the charge of a
# single electron moving across an electric potential difference of one volt.
electron_volt	= 1.60217662 * 10**(-19)	# joules kg*m^2/s^2

# Length dictionary
length_dictionary 	= {"fermi" : 10**(-15), "nano" : 10**(-9), "micro" : 10**(-6)}

# Energy dictionary
energy_dictionary 	= {"mev" : 10**6, "gev" : 10**9}

# LENGTHS IN EV: hbarc in different unit lengths
hbarc_meters_ev		= hbarc / electron_volt 						# [eV] meters
hbarc_micro_ev		= hbarc_meters_ev / length_dictionary["micro"]	# [MeV] micro-meters
hbarc_nano_ev		= hbarc_meters_ev / length_dictionary["nano"]	# [GeV] nano-meters
hbarc_fermi_ev		= hbarc_meters_ev / length_dictionary["fermi"]	# [PeV] fermi-meters

# ENERGY: hbarc in different energy units
hbarc_meters_mev 	= hbarc_meters_ev / energy_dictionary["mev"]	# MeV micro
hbarc_meters_gev 	= hbarc_meters_ev / energy_dictionary["gev"]	# GeV meters

hbarc_micro_mev		= hbarc_meters_mev / length_dictionary["micro"]	# MeV micro
hbarc_micro_gev		= hbarc_meters_gev / length_dictionary["micro"]	# GeV micro

hbarc_nano_mev		= hbarc_meters_mev / length_dictionary["nano"] 	# MeV nano
hbarc_nano_gev		= hbarc_meters_gev / length_dictionary["nano"]  # GeV nano

hbarc_fermi_mev		= hbarc_meters_mev / length_dictionary["fermi"] # MeV fermi
hbarc_fermi_gev		= hbarc_meters_gev / length_dictionary["fermi"]	# GeV fermi

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

def convert_meters_to_ev(fm):
	return hbarc*fm

def convert_ev_to_meters(fm):
	return fm

print getLatticeSpacing(6.45) * length_dictionary["fermi"] * 48
print "Lattice spacing: %.4f [fermi]" % getLatticeSpacing(6.0)
print "Energy scale:	%.4f [GeV]" % (getLatticeSpacing(6.45) / hbarc_fermi_gev)
print hbarc_nano_mev
print hbarc_fermi_gev
print hbarc_micro_ev

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
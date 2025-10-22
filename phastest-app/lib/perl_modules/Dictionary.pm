#!/usr/bin/perl
# April 9, 2011 Rah. This library is a filter to pick up phage like proteins by names shown in
# annoatations (in genbank files or blast results)

# dictonary

my @general_key_word = (
	"phage", "prophage", "bacteriophage", 
	"lambda", "T1", "T3", "T4", "T5", "T7", "Mu ", 
	"procapsid", "capsid", "endolysin", "lysis", "tail", "lytic",
	"host", "fibre", "holin", "scaffold", "lysozyme", "integrase", "terminase",
	"injection", "sheath", "virulence", "viral", "fimbrillin", "packag", "head",
	"portal", "excisionase", "coat", "toxin");

my @strict_key_word = (
	"phage", "prophage", "bacteriophage", 
	"lambda", "T1", "T3", "T4", "T5", "T7", "Mu ", 
	"procapsid", "capsid", "endolysin", "lysis", "fibre", "holin", "coat", "head",
	"host",  "lysozyme", "integrase", "terminase", "excisionase", "packaging",
	"injection", "sheath", "virulence", "fimbrillin", "portal");

my @strict_exclude_word = ("transposase", "inhibit");
	
my @phage_names =  ("APSE-1", "BIP-1", "P1201", "BcepF1", "VfO3K6", "phiadh", "EPS7", "RB43"
	, "PRD1", "Spud", "KS9", "PRR1", "phi13", "BK5-T", "PVL", "Gumball"
	, "IEBH", "L2", "Acj9", "phiETA2", "K139", "MP38", "52A", "phiYS40"
	, "Orion", "PA11", "phiEcoM", "PY54", "phi8", "kappa", "AP50", "MV-L1"
	, "PH10", "His1", "O1205", "Ac42", "VP93", "PM2", "PLot", "phiP27"
	, "Boomer", "phiL7", "KP32", "PG1", "Lv-1", "LUZ7", "KSF-1phi", "VT2-Sakai"
	, "bIL312", "MP29", "DT1", "Cd", "2638A", "Ike", "Pukovnik", "P23-77"
	, "Solon", "SIO1", "CrimD", "M102", "phig1e", "S1", "phiKZ", "MM1"
	, "S-PM2", "Pipefish", "Angel", "Tuc2009", "HK022", "tp310-3", "Yepe2", "KC5a"
	, "F108", "Lc-Nu", "P4", "42E", "phiPVL-CN125", "SPO1", "YYZ-2008", "SP6"
	, "L-413C", "BP-4795", "APSE-2", "phiKO2", "K1-5", "N15", "TLS", "ET08"
	, "S-RSM4", "RM378", "phi105", "ST64T", "EJ-1", "A118", "44RR2.8t", "A2"
	, "Min27", "Bxz1", "BMP-1", "BcepC6B", "phiEf11", "JK06", "Aaphi23", "fs2"
	, "M13", "Halo", "phiMR25", "16-3", "MmP1", "D108", "P-SSM4", "CNPH82"
	, "phi-2", "A511", "vB_EcoM-VR7", "Bcep1", "phiE12-2", "HK97", "P68"
	, "M6", "GIL16c", "phiYeO3-12", "phiSLT", "tp310-2", "phiJL-1", "Phlyer"
	, "EE36phi1", "ID2", "phiO18P", "PT2", "bIL310", "phiEa21-4", "Lb338-1", "Felix"
	, "JS10", "VfO4K68", "Qbeta", "DSS3phi2", "SS3e", "Xop411", "JS98"
	, "VGJphi", "phiLC3", "PaP3", "Acj61", "TP901-1", "14-1", "KSY1"
	, "EL", "phBC6A51", "P1", "K11", "cdtI", "ES18", "r1t", "phiSMA9"
	, "HP2", "KP34", "Bcep781", "P087", "phiHSIC", "phi52237", "RB51", "Bcep176"
	, "Giles", "phiPLPE", "B3", "Kamchatka", "T5", "jj50", "phiC2", "39-O"
	, "1-R8A2B", "Xfas53", "EW", "tp310-1", "Bxz2", "Brujita", "3A", "RSL1"
	, "B054", "Porky", "BA3", "ALQ13.2", "F10", "VP2", "DD5", "44AHJD"
	, "PAJU2", "Sf6", "LUZ24", "phiNM3", "RSM1", "SkV1_CR2-3x", "P2"
	, "Berlin", "T4", "RSM3", "P9", "P22", "Lockley", "Ragged Hills", "ul36"
	, "PB1", "YuA", "LeBron", "E2", "TP21-L", "Xp10", "STSV1", "Phieco32"
	, "epsilon34", "phiFL3A", "KBG", "Lrm1", "phiAT3", "PBC5", "Angelica", "SETP3"
	, "phi29", "EFAP-1", "WA13 sensu lato", "phikF77", "G1", "SP18", "Butterscotch", "c2"
	, "Sfi21", "Pf-WMP3", "VP882", "Chah", "phiPVL108", "Bam35c", "PBI1", "phiEF24C"
	, "BZ13", "bIL311", "Catera", "phiA1122", "MP22", "phiMH2K", "L5"
	, "Pf3", "11b", "bIBB29", "KP15", "phi2958PVL", "phi1026b", "Min1", "JSE"
	, "SO-1", "D29", "BJ1", "SVTS2", "phiSM101", "933W", "CPAR39", "Q54"
	, "Kostya", "bIL170", "ROSA", "phiHAP-1", "VSK", "P008", "K1E", "asccphi28"
	, "phiAS4", "Cp-1", "Aeh1", "CC31", "Vi II-E1", "VP5", "T3", "ST64B"
	, "P35", "PA6", "119X", "Lj965", "Phaedrus", "Myrna", "CMP1", "Wildcat"
	, "Chp1", "phiRSA1", "Bcep22", "phiE255", "phi CD119", "B103", "phiMHaA1", "HF2"
	, "N4", "B025", "Rosebush", "F116", "WV8", "Lj771", "BA14", "Fels-1"
	, "IME08", "Bethlehem", "P335", "bIL309", "phi3626", "BFK20", "I2-2", "Pf-WMP4"
	, "KVP40", "G4 sensu lato", "Che12", "P40", "Abc2", "phiSASD1", "SE1", "A006"
	, "If1", "Chp2", "Ma-LMM01", "ID18 sensu lato", "B5", "Gifsy-2", "Che9d", "phiCTX"
	, "phiAsp2", "phiFL2A", "bIL286", "SN", "bIL67", "T7", "psiM2", "phi2954"
	, "psiM100", "Tweety", "K1F", "phiCTP1", "Cooper", "alpha3", "BPP-1", "BcepB1A"
	, "Konstantine", "BcepIL02", "bIL285", "phiKMV", "SPbeta", "SPP1", "EcoDS1", "PhiCh1"
	, "HF1", "IN93", "syn9", "BcepGomr", "0305phi8-36", "Jasper", "PT1028", "FI sensu lato"
	, "Pacc40", "Fah", "phiMFV1", "RSS1", "Twort", "GBSV1", "ST104", "Peaches"
	, "Pf1", "RB32", "phiX174", "PH15", "c-st", "LP65", "1-C74", "Rizal"
	, "P74-26", "Era103", "PP7", "phi 12", "P-SSP7", "LKA1", "BcepNazgul"
	, "phiAS5", "PMC", "Gifsy-1", "SSL-2009a", "RTP", "Adjutor", "Troll4", "phiSauS-IPLA88"
	, "RSB1", "201phi2-1", "80alpha", "HP1", "Llij", "gh-1", "Cf1c", "Vf12"
	, "phi6", "VWB", "PT5", "C1", "Ramsey", "SAP-2", "Lj928", "RB14"
	, "BPs", "LMA2", "Sulfolobus turreted icosahedral", "phiN315", "T1", "AP205", "BCJA1c", "phiSboM-AG3"
	, "fs1", "Omega", "p12J", "SfV", "SMP", "P23-45", "SAP-26", "13a"
	, "phiSauS-IPLA35", "LBL3", "P-SSM2", "Av-1", "Predator", "phiETA3", "OP2", "S1249"
	, "Vf33", "phiW-14", "SM1", "Fruitloop", "His2", "DMS3", "VHML", "phiCPG1"
	, "phiBT1", "Mx8", "KS10", "rv5", "OP1", "Bcep43", "RB69", "TM4"
	, "phiC31", "LUZ19", "PSS2", "phBC6A52", "phiMR11", "GA-1", "phi3396", "PsP3"
	, "Sfi19", "phiV10", "Phi1", "Ardmore", "Corndog", "Nigel", "A500", "sk1"
	, "LKD16", "D3112", "F8", "RB16", "phiE202", "phiFL4A", "phi644-2", "VpV262"
	, "B40-8", "Fels-2", "Cjw1", "lambda", "RB49", "X2", "Xp15", "phiPV83"
	, "phi12", "MAV1", "U2", "Gamma", "Qyrzula", "WBeta", "Che8", "phiETA"
	, "VP4", "Mu", "Syn5", "BcepNY3", "epsilon15", "c341", "MS2", "phiSG-JL2"
	, "P954", "phiSG1", "VEJphi", "phiNIH1.1", "PaP2", "phiCD27", "Barnyard", "Cali"
	, "Sfi11", "D3", "St-1", "LIT1", "BcepMu", "Che9c", "phiNM", "phiJL001"
	, "phiE125", "HK620", "Kvp1", "phiFL1A", "P60", "Bxb1", "LL-H"
	, "mu1/6", "ScottMcG"
);
# this function 
sub generalRelated {
	my $name = shift;
	foreach (@general_key_word){
		if ($name =~m/\b$_/i){
			return 1;
		}
	}
#	foreach (@phage_names){
#		if ($name =~m/\b$_/i){
#			return 1;
#		}
#	}
	return 0;
}

sub strictRelated {
	my $name = shift;
	my $ok = 0;
	foreach (@strict_key_word){
		if ($name =~m/\b$_/i){
			$ok = 1;
			last;
		}
	}
	if ($ok == 1){
		foreach (@strict_exclude_word){
			if ($name =~m/\b$_/i){
				return 0;
			}
		}
		return 1;
	}
	return 0;
}



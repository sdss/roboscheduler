#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <stdio.h>
#include <vector>
#include <pybind11/pybind11.h> // must be first
#include <string>
#include "cadenceCore.h"

namespace py = pybind11;

CadenceCore::CadenceCore(std::string name,
												 int nepochs,
												 Instrument instrument,
												 std::vector<float> skybrightness,
												 std::vector<float> delta,
												 std::vector<float> delta_min,
												 std::vector<float> delta_max,
												 std::vector<int> nexp) :
	name(name), nepochs(nepochs), instrument(instrument),
	skybrightness(skybrightness), delta(delta), delta_min(delta_min),
	delta_max(delta_max), nexp(nexp)
{
	
	// Count up the exposures for convenience, set up "epoch_indx"
	// from epoch number into exposure list (note it has nepochs+1 entries),
	// and set up "epochs" giving epoch number for each exposure
	nexp_total = 0;
	for(auto i = 0; i < nepochs; i++) {
		epoch_indx.push_back(nexp_total);
		nexp_total += nexp[i];
		for(auto j = 0; j < nexp[i]; j++)
			epochs.push_back(i);
	}
	epoch_indx.push_back(nexp_total);
}

std::string CadenceCore::__str__()
{
	return(epochText());
}

std::string CadenceCore::epochText()
{
	std::string out, tmps;
	char tmp[2000];

	out = name + "\n";
	out = out + " nepochs=" + std::to_string(nepochs) + "\n";

	if(instrument == ApogeeInstrument)
		out = out + " instrument=APOGEE\n";
	else
		out = out + " instrument=BOSS\n";

	out = out + " skybrightness=";
	for(auto i = 0; i < nepochs; i++) {
		sprintf(tmp, " %4.2f", skybrightness[i]);
		out = out + tmps.assign(tmp);
	}
	out = out + "\n";

	out = out + " delta=";
	for(auto i = 0; i < nepochs; i++) {
		sprintf(tmp, " %4.2f", delta[i]);
		out = out + tmps.assign(tmp);
	}
	out = out + "\n";

	out = out + " delta_min=";
	for(auto i = 0; i < nepochs; i++) {
		sprintf(tmp, " %4.2f", delta_min[i]);
		out = out + tmps.assign(tmp);
	}
	out = out + "\n";

	out = out + " delta_max=";
	for(auto i = 0; i < nepochs; i++) {
		sprintf(tmp, " %4.2f", delta_max[i]);
		out = out + tmps.assign(tmp);
	}
	out = out + "\n";

	out = out + " nexp=";
	for(auto i = 0; i < nepochs; i++) {
		sprintf(tmp, " %d", nexp[i]);
		out = out + tmps.assign(tmp);
	}

	return(out);
}

bool CadenceCore::epochsConsistency(CadenceCore target_cadence,
																		std::vector<int> epochs) {
	int ok;
	float dtotmin, dtotmax;

	ok = (target_cadence.nexp[0] <= nexp[epochs[0]]) &
		(target_cadence.skybrightness[0] >= skybrightness[epochs[0]]);
	if(!ok)
		return(false);

	for(unsigned long i = 1; i < epochs.size(); i++) {
		if(target_cadence.delta[i] >= 0) {
			dtotmin = 0.;
			for(auto j = epochs[i - 1] + 1; j <= epochs[i]; j++)
				dtotmin += delta_min[j];
			dtotmax = 0.;
			for(auto j = epochs[i - 1] + 1; j <= epochs[i]; j++)
				dtotmax += delta_max[j];
			ok = (target_cadence.delta_min[i] <= dtotmin) &
				(target_cadence.delta_max[i] >= dtotmax) &
				(target_cadence.nexp[i] <= nexp[epochs[i]]) &
				(target_cadence.skybrightness[i] >= skybrightness[epochs[i]]);
		} else {
			ok = (target_cadence.nexp[i] <= nexp[epochs[i]]) &
				(target_cadence.skybrightness[i] >= skybrightness[epochs[i]]);
		}
		if(!ok)
			return(false);
	}
	return(true);
}

std::vector<std::vector<int>> CadenceCore::cadenceConsistency(CadenceCore target_cadence) {
	std::vector<std::vector<int>> epochs_list;
	std::vector<std::vector<int>> current_epochs_list;
	std::vector<int> epochs;
	int ok, ifield_start;

	// Initialize solutions to those to start with
	for(auto istart = 0; istart < nepochs; istart++) {
		epochs.clear();
		epochs.push_back(istart);
		if(epochsConsistency(target_cadence, epochs))
			epochs_list.push_back(epochs);
	}

	// Return if nowhere to begin
	if(epochs_list.size() == 0) {
		return(epochs_list);
	}

	// Successively find solutions through each target epoch
	for(int itarget = 1; itarget < target_cadence.nepochs; itarget++) {
		current_epochs_list = epochs_list;
		epochs_list.clear();

		// Each time, check all previous working solutions
		for(unsigned long i = 0; i < current_epochs_list.size(); i++) {
			epochs = current_epochs_list[i];
			ifield_start = epochs.back() + 1;

			// And find all subsequent field epochs that will work for each
			for(int ifield = ifield_start; ifield < nepochs; ifield++) {
				epochs = current_epochs_list[i];
				epochs.push_back(ifield);
				ok = epochsConsistency(target_cadence, epochs);
				if(ok) {
					epochs_list.push_back(epochs);
				}
			}
			epochs.clear();
		}

	}

	return(epochs_list);
	
}

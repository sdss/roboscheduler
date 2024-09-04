#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <stdio.h>
#include <vector>
#include <pybind11/pybind11.h> // must be first
#include <pybind11/numpy.h> // must be first
#include <string>
#include "cadenceCore.h"

namespace py = pybind11;

CadenceCore::CadenceCore(std::string name,
												 int nepochs,
												 py::array_t<float> skybrightness,
												 py::array_t<float> delta,
												 py::array_t<float> delta_min,
												 py::array_t<float> delta_max,
												 py::array_t<int> nexp,
												 py::array_t<float> max_length,
												 py::array_t<float> min_moon_sep,
												 py::array_t<float> min_deltav_ks91,
												 py::array_t<float> min_twilight_ang,
												 py::array_t<float> max_airmass,
												 py::array_t<int> epoch_indx,
												 py::array_t<int> epochs) :
	name(name), nepochs(nepochs),
	skybrightness(skybrightness), delta(delta), delta_min(delta_min),
	delta_max(delta_max), nexp(nexp), max_length(max_length),
	min_moon_sep(min_moon_sep), min_deltav_ks91(min_deltav_ks91),
	min_twilight_ang(min_twilight_ang), max_airmass(max_airmass),
	epoch_indx(epoch_indx), epochs(epochs)
{
	int *nexp_a = (int *) nexp.request().ptr;
	int *epoch_indx_a = (int *) epoch_indx.request().ptr;
	int *epochs_a = (int *) epochs.request().ptr;

	// Count up the exposures for convenience, set up "epoch_indx"
	// from epoch number into exposure list (note it has nepochs+1 entries),
	// and set up "epochs" giving epoch number for each exposure
	nexp_total = 0;
	for(auto i = 0; i < nepochs; i++) {
		epoch_indx_a[i] = nexp_total;
		for(auto j = 0; j < nexp_a[i]; j++)
			epochs_a[nexp_total + j] = i;
		nexp_total += nexp_a[i];
	}
	epoch_indx_a[nepochs] = nexp_total;
}

std::string CadenceCore::__str__()
{
	return(epochText());
}

std::string CadenceCore::epochText()
{
	std::string out, tmps;
	char tmp[2000];

	int *nexp_a = (int *) nexp.request().ptr;
	float *max_length_a = (float *) max_length.request().ptr;
	float *min_moon_sep_a = (float *) min_moon_sep.request().ptr;
	float *min_deltav_ks91_a = (float *) min_deltav_ks91.request().ptr;
	float *min_twilight_ang_a = (float *) min_twilight_ang.request().ptr;
	float *max_airmass_a = (float *) max_airmass.request().ptr;
	float *skybrightness_a = (float *) skybrightness.request().ptr;
	float *delta_a = (float *) delta.request().ptr;
	float *delta_min_a = (float *) delta_min.request().ptr;
	float *delta_max_a = (float *) delta_max.request().ptr;

	out = "[" + name + "]\n";
	out = out + " nepochs=" + std::to_string(nepochs) + "\n";

	out = out + " skybrightness=";
	for(auto i = 0; i < nepochs; i++) {
		sprintf(tmp, " %4.2f", skybrightness_a[i]);
		out = out + tmps.assign(tmp);
	}
	out = out + "\n";

	out = out + " delta=";
	for(auto i = 0; i < nepochs; i++) {
		sprintf(tmp, " %4.2f", delta_a[i]);
		out = out + tmps.assign(tmp);
	}
	out = out + "\n";

	out = out + " delta_min=";
	for(auto i = 0; i < nepochs; i++) {
		sprintf(tmp, " %4.2f", delta_min_a[i]);
		out = out + tmps.assign(tmp);
	}
	out = out + "\n";

	out = out + " delta_max=";
	for(auto i = 0; i < nepochs; i++) {
		sprintf(tmp, " %4.2f", delta_max_a[i]);
		out = out + tmps.assign(tmp);
	}
	out = out + "\n";

	out = out + " nexp=";
	for(auto i = 0; i < nepochs; i++) {
		sprintf(tmp, " %d", nexp_a[i]);
		out = out + tmps.assign(tmp);
	}
	out = out + "\n";

	out = out + " max_length=";
	for(auto i = 0; i < nepochs; i++) {
		sprintf(tmp, " %4.2f", max_length_a[i]);
		out = out + tmps.assign(tmp);
	}
	out = out + "\n";

	out = out + " min_moon_sep=";
	for(auto i = 0; i < nepochs; i++) {
		sprintf(tmp, " %4.2f", min_moon_sep_a[i]);
		out = out + tmps.assign(tmp);
	}
	out = out + "\n";

	out = out + " min_deltav_ks91=";
	for(auto i = 0; i < nepochs; i++) {
		sprintf(tmp, " %4.2f", min_deltav_ks91_a[i]);
		out = out + tmps.assign(tmp);
	}
	out = out + "\n";

	out = out + " min_twilight_ang=";
	for(auto i = 0; i < nepochs; i++) {
		sprintf(tmp, " %4.2f", min_twilight_ang_a[i]);
		out = out + tmps.assign(tmp);
	}
	out = out + "\n";

	out = out + " max_airmass=";
	for(auto i = 0; i < nepochs; i++) {
		sprintf(tmp, " %4.2f", max_airmass_a[i]);
		out = out + tmps.assign(tmp);
	}
	out = out + "\n";

	out = out + "\n";

	return(out);
}

bool CadenceCore::epochsConsistency(CadenceCore target_cadence,
																		std::vector<int> epochs,
																		bool skybrightnessOnly,
																		bool inorder) {

	int ok;
	std::vector<int> target_epochs;

	// If we are supposed to sample the target_cadence the same way
	// as the object cadence.
	if(inorder) {
		for(unsigned long i = 0; i < epochs.size(); i++) {
			target_epochs.push_back(epochs[i]);
		} // end for 
	} else {
		for(unsigned long i = 0; i < epochs.size(); i++) {
			target_epochs.push_back(i);
		} // end for 
	} // end if..else

	ok = specificEpochsConsistency(target_cadence, epochs, target_epochs, skybrightnessOnly);

	return(ok);
}

bool CadenceCore::specificEpochsConsistency(CadenceCore target_cadence,
																						std::vector<int> epochs,
																						std::vector<int> target_epochs,
																						bool skybrightnessOnly) {

	int ok;
	float dtotmin, dtotmax;
  std::vector<int> nexp_count;

	nexp_count.resize(nepochs);
	for(int i = 0; i < nepochs; i++)
		nexp_count[i] = 0;

	int *t_nexp_a = (int *) target_cadence.nexp.request().ptr;
	int *nexp_a = (int *) nexp.request().ptr;

	float *t_skybrightness_a = (float *) target_cadence.skybrightness.request().ptr;
	float *skybrightness_a = (float *) skybrightness.request().ptr;

	float *t_min_deltav_ks91_a = (float *) target_cadence.min_deltav_ks91.request().ptr;
	float *t_max_airmass_a = (float *) target_cadence.max_airmass.request().ptr;
	float *t_delta_a = (float *) target_cadence.delta.request().ptr;
	float *t_delta_min_a = (float *) target_cadence.delta_min.request().ptr;
	float *t_delta_max_a = (float *) target_cadence.delta_max.request().ptr;
	float *min_deltav_ks91_a = (float *) min_deltav_ks91.request().ptr;
	float *max_airmass_a = (float *) max_airmass.request().ptr;
	float *delta_min_a = (float *) delta_min.request().ptr;
	float *delta_max_a = (float *) delta_max.request().ptr;
	float *delta_a = (float *) delta.request().ptr;

	if(skybrightnessOnly)
		ok = (t_skybrightness_a[target_epochs[0]] >= skybrightness_a[epochs[0]]);
	else
		ok = (t_min_deltav_ks91_a[target_epochs[0]] <= min_deltav_ks91_a[epochs[0]]) &
			(t_max_airmass_a[target_epochs[0]] >= max_airmass_a[epochs[0]]);

	ok = ok & (t_nexp_a[target_epochs[0]] <= nexp_a[epochs[0]] - nexp_count[epochs[0]]);

	if(!ok)
		return(false);

	nexp_count[epochs[0]] += t_nexp_a[target_epochs[0]];

	for(unsigned long i = 1; i < epochs.size(); i++) {
		if(skybrightnessOnly)
			ok = (t_skybrightness_a[target_epochs[i]] >= skybrightness_a[epochs[i]]);
		else
			ok = (t_min_deltav_ks91_a[target_epochs[i]] <= min_deltav_ks91_a[epochs[i]]) &
				(t_max_airmass_a[target_epochs[i]] >= max_airmass_a[epochs[i]]);
		if(!ok)
			return(false);

		if(t_delta_a[target_epochs[i]] > 0) {
			dtotmin = 0.;
			for(auto j = epochs[i - 1] + 1; j <= epochs[i]; j++)
				dtotmin += delta_min_a[j];
			dtotmax = 0.;
			for(auto j = epochs[i - 1] + 1; j <= epochs[i]; j++)
				dtotmax += delta_max_a[j];
			ok = ok & (t_delta_min_a[target_epochs[i]] <= dtotmin) &
				(t_delta_max_a[target_epochs[i]] >= dtotmax) &
				(t_nexp_a[target_epochs[i]] <= nexp_a[epochs[i]] - nexp_count[epochs[i]]);
		} else {
			ok = ok & (t_nexp_a[target_epochs[i]] <= nexp_a[epochs[i]] - nexp_count[epochs[i]]);
		}
		if(!ok)
			return(false);

		nexp_count[epochs[i]] += t_nexp_a[target_epochs[i]];
	}
	return(true);
}


std::vector<std::vector<int>> CadenceCore::cadenceConsistency(CadenceCore target_cadence,
																															bool skybrightnessOnly) {
	std::vector<std::vector<int>> epochs_list;
	std::vector<int> target_epochs;

	for(auto i = 0; i < target_cadence.nepochs; i++)
		target_epochs.push_back(i);

	epochs_list = specificCadenceConsistency(target_cadence, target_epochs,
																					 skybrightnessOnly, 0, false);

	return(epochs_list);
}


std::vector<std::vector<int>> CadenceCore::specificCadenceConsistency(CadenceCore target_cadence,
																																			std::vector<int> target_epochs,
																																			bool skybrightnessOnly,
																																			long unsigned int limit,
																																			bool sequential) {
	std::vector<std::vector<int>> epochs_list;
	std::vector<std::vector<int>> current_epochs_list;
	std::vector<int> current_target_epochs;
	std::vector<int> current_epochs;
	int ok, ifield_start, ifield_end, nexp_left_target, nexp_left_field, target_nepochs;

	target_nepochs = target_epochs.size();

	int *t_nexp_a = (int *) target_cadence.nexp.request().ptr;
	int *nexp_a = (int *) nexp.request().ptr;

	// Initialize solutions to those to start with
	for(auto istart = 0; istart < nepochs; istart++) {
		current_epochs.clear();
		current_epochs.push_back(istart);
		if(specificEpochsConsistency(target_cadence, current_epochs, target_epochs,
																 skybrightnessOnly)) {
			// Check number of exposures left total; if the field
			// doesn't have enough exposures left, no point in continuing
			// the check
			nexp_left_target = 0;
			for(auto j = 1; j < target_nepochs; j++)
				nexp_left_target += t_nexp_a[target_epochs[j]];
			nexp_left_field = nexp_a[istart] - t_nexp_a[target_epochs[0]];
			for(auto j = istart + 1; j < nepochs; j++)
				nexp_left_field += nexp_a[j];

			if(nexp_left_field >= nexp_left_target) {
				epochs_list.push_back(current_epochs);
			}
		}
	}

	// Return if nowhere to begin
	if(epochs_list.size() == 0) {
		return(epochs_list);
	}

	// Successively find solutions through each target epoch
	for(int indx = 1; indx < target_nepochs; indx++) {
		current_epochs_list = epochs_list;
		epochs_list.clear();

		// Each time, check all previous working solutions
		for(unsigned long i = 0; i < current_epochs_list.size(); i++) {
			current_epochs = current_epochs_list[i];
			ifield_start = current_epochs.back(); // allow to repeat the epoch

			// Just check the next one if the input is sequential (and sequential option on)
			if((sequential) & (target_epochs[indx] <= target_epochs[indx - 1] + 1)) {
				ifield_end = ifield_start + 1;
				if(ifield_end > nepochs -1)
					ifield_end = nepochs - 1;
			} else {
				ifield_end = nepochs - 1;
			}

			// Find subsequent field epochs that will work for each
			for(int ifield = ifield_start; ifield <= ifield_end; ifield++) {
				current_epochs = current_epochs_list[i];
				current_epochs.push_back(ifield);
				ok = specificEpochsConsistency(target_cadence, current_epochs, target_epochs,
																			 skybrightnessOnly);
				if(ok) {
					// Check number of exposures left total; if the field
					// doesn't have enough exposures left, no point in
					// continuing the check; note that this check includes any
					// previously used exposures of the current epoch, even though
					// some of those may have been used already. That will allow
					// a non-viable case through, but the epochConsistency()
					// check will catch those issues in subsequent target epochs.
					nexp_left_target = 0;
					for(int j = current_epochs.size(); j < target_nepochs; j++)
						nexp_left_target += t_nexp_a[target_epochs[j]];
					nexp_left_field = nexp_a[ifield] - t_nexp_a[target_epochs[current_epochs.size() - 1]];
					for(auto j = ifield + 1; j < nepochs; j++)
						nexp_left_field += nexp_a[j];
					if(nexp_left_field >= nexp_left_target)
						epochs_list.push_back(current_epochs);
				}
			}
			current_epochs.clear();

			// If we have set a limit, don't try to get more than that. Note that
			// this limits how many options there are at any stage, so in the end
			// there might be fewer than the limit if some of the options do not
			// pan out.
			if((limit > 0) && (epochs_list.size() > limit))
				break;
		}

	}

	return(epochs_list);
	
}

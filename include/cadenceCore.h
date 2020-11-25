# pragma once
#include <vector>
#include <map>

enum Instrument {ApogeeInstrument, BossInstrument};

class CadenceCore {
public:
	CadenceCore(std::string name, int nepochs, Instrument instrument,
							std::vector<float> skybrightness,
							std::vector<float> delta,
							std::vector<float> delta_min, std::vector<float> delta_max,
							std::vector<int> nexp);
	std::string __str__();
	std::string epochText();
	bool epochsConsistency(CadenceCore target_cadence,
												 std::vector<int> epochs);
	std::vector<std::vector<int>> cadenceConsistency(CadenceCore target_cadence);
	std::string name;
	int nepochs;
	int nexp_total;
	Instrument instrument;
	std::vector<float> skybrightness;
	std::vector<float> delta;
	std::vector<float> delta_min;
	std::vector<float> delta_max;
	std::vector<int> nexp;
	std::vector<int> epoch_indx;
	std::vector<int> epochs;
};

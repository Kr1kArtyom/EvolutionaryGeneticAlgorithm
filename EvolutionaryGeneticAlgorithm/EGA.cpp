#include "EGA.h"
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <map>

Individ::Individ(int numCities, float _probabilityGeneMutation) : encoding(numCities), adaptability(0), probabilityGeneMutation(_probabilityGeneMutation) {
	for (int i = 0; i < numCities; i++) {
		encoding[i] = i;
	}
}

float Individ::CalculateAdaptability(const MatrixDistances distances) {
	adaptability = 0.00;
	for (int i = 1; i < encoding.size(); i++) {
		int src = encoding[i - 1], dest = encoding[i];
		adaptability += distances[src][dest];
	}
	adaptability = (float)distances.size() * 15 - adaptability;
	return adaptability;
}

void Individ::SwapGenes(int indFirst, int indSecond) {
	int temp = this->encoding[indFirst];
	this->encoding[indFirst] = this->encoding[indSecond];
	this->encoding[indSecond] = temp;
}

Roulette::Roulette(std::vector<float> adaptabilities) {
	sizeRoulette = 0;
	for (int i = 0; i < adaptabilities.size(); i++) {
		sizeRoulette += adaptabilities[i];
		roulette.push_back(sizeRoulette);
	}
}

int Roulette::ChooseObjectByRoulette() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> randInt(0, sizeRoulette);

	int randIndivid = randInt(gen);
	for (int i = 0; i < roulette.size(); i++) {
		if (randIndivid <= roulette[i]) return i;
	}
}

EGA::EGA(int _sizePopulation, MatrixDistances _distances, int _maxNumIterations, int _maxNumIterationsWithoutImprovements) : sizePopulation(_sizePopulation),
numGenes(_distances.size()), distances(_distances), population(_sizePopulation), maxNumIterations(_maxNumIterations), numIterations(0), BestIndivid(0),
maxNumIterationsWithoutImprovements(_maxNumIterationsWithoutImprovements), NumIterationsWithoutImprovements(0) {}

void EGA::swapIndivids(Individ& indFirst, Individ& indSecond) {
	Individ temp = indFirst;
	indFirst = indSecond;
	indSecond = temp;
}

Population EGA::generateRandomPopulations() {
	for (int it = 0; it < sizePopulation; it++) {
		Individ newIndivid(numGenes);
		std::random_device rd;
		std::mt19937 g(rd());
		shuffle(newIndivid.encoding.begin(), newIndivid.encoding.end(), g);
		newIndivid.CalculateAdaptability(distances);
		population[it] = newIndivid;
	}
	return population;
}

Population EGA::generateInitialPopulation() {
	std::random_device rd;
	std::mt19937 gen(rd());

	for (int it = 0; it < sizePopulation; it++) {
		Individ newIndivid(distances.size());
		int curIndCity = 0;
		int nextCity = std::uniform_int_distribution<int>(0, distances.size() - 1)(gen);
		newIndivid.SwapGenes(newIndivid.encoding[curIndCity], nextCity);

		while (curIndCity != numGenes) {
			auto nextCity = std::min_element(newIndivid.encoding.begin() + curIndCity + 1, newIndivid.encoding.end(), [&](int i1, int i2) {
				return distances[newIndivid.encoding[curIndCity]][i1] < distances[newIndivid.encoding[curIndCity]][i2];
				});
			newIndivid.encoding[++curIndCity] = *nextCity;
		}

		newIndivid.CalculateAdaptability(distances);
		population.push_back(newIndivid);
	}
	return population;
}

Population EGA::generateGreedyPopulations() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> randInt(0, numGenes - 1);
	for (int it = 0; it < sizePopulation; it++) {
		Individ newIndivid(numGenes);
		int curNumCities = 0;
		int nextCity = randInt(gen);
		newIndivid.SwapGenes(newIndivid.encoding[curNumCities], nextCity);
		curNumCities++;

		while (curNumCities != numGenes) {
			std::vector<std::pair<int, int>> Z;
			for (auto iter = newIndivid.encoding.begin(); iter != newIndivid.encoding.begin() + curNumCities; iter++) {
				auto it = std::min_element(newIndivid.encoding.begin() + curNumCities, newIndivid.encoding.end(), [&](int i1, int i2) {
					return distances[*iter][i1] < distances[*iter][i2];
					});
				Z.push_back(std::make_pair(*iter, *it));
			}

			auto z = min_element(Z.begin(), Z.end(), [&](std::pair<int, int> i1, std::pair<int, int> i2) {
				return distances[i1.first][i1.second] < distances[i2.first][i2.second];
				});

			nextCity = (*z).second;
			newIndivid.SwapGenes(newIndivid.encoding[curNumCities], nextCity);
			curNumCities++;
		}
		newIndivid.CalculateAdaptability(distances);
		population[it] = newIndivid;
	}
	return population;
}

Population EGA::RandomSelection() {
	std::random_device rd;
	std::mt19937 g(rd());
	shuffle(population.begin(), population.end(), g);
	return population;
}

Population EGA::PositiveAssortativeMating() {
	std::sort(population.begin(), population.end(), [](const Individ& a, const Individ& b) {
		return a.adaptability > b.adaptability;
		});
	return population;
}

Population EGA::NegativeAssortativeMating() {
	std::sort(population.begin(), population.end(), [](const Individ& a, const Individ& b) {
		return std::abs(a.adaptability - b.adaptability) > std::abs(a.adaptability - b.adaptability);
		});
	return population;
}

Population EGA::InbreedingSelection() {
	for (int i = 0; i < sizePopulation / 2 * 2; i += 2) {
		int curIndivid = i;
		int maxMatches = 0;
		int indexMaxMatches = 0;
		for (int j = i + 1; j < sizePopulation; j++) {
			int matches = 0;
			for (int k = 0; k < numGenes; k++)
				if (population[curIndivid].encoding[k] == population[j].encoding[k]) matches++;
			if (matches > maxMatches) {
				maxMatches = matches;
				indexMaxMatches = j;
			}
		}
		swapIndivids(population[curIndivid + 1], population[indexMaxMatches]);
	}
	return population;
}

Population EGA::OutbreedingSelection() {
	for (int i = 0; i < sizePopulation / 2 * 2; i += 2) {
		int curIndivid = i;
		int minMatches = numGenes;
		int indexMinMatches = 0;
		for (int j = i + 1; j < sizePopulation; j++) {
			int matches = 0;
			for (int k = 0; k < numGenes; k++)
				if (population[curIndivid].encoding[k] == population[j].encoding[k]) matches++;
			if (matches < minMatches) {
				minMatches = matches;
				indexMinMatches = j;
			}
		}
		swapIndivids(population[curIndivid + 1], population[indexMinMatches]);
	}
	return population;
}

Population EGA::OrderCrossover() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> dist(0, numGenes - 1);

	for (int it = 0; it < sizePopulation / 2 * 2; it++) {
		Individ descendant(numGenes);
		Individ parentFir = population[it];
		Individ parentSec;
		if (it % 2 != 0) parentSec = population[it - 1];
		else parentSec = population[it + 1];

		int startPos = dist(gen);
		int endPos = dist(gen);
		if (endPos < startPos) {
			std::swap(startPos, endPos);
		}

		std::unordered_set<int> selected;
		for (int i = startPos; i <= endPos; ++i) {
			descendant.encoding[i] = parentFir.encoding[i];
			selected.insert(parentFir.encoding[i]);
		}

		int j = 0;
		for (int i = 0; i < numGenes; ++i) {
			if (i < startPos || i > endPos) {
				while (selected.find(parentSec.encoding[j]) != selected.end()) {
					++j;
				}
				descendant.encoding[i] = parentSec.encoding[j];
				++j;
			}
		}
		descendant.CalculateAdaptability(distances);
		descendants.push_back(descendant);
	}
	if (sizePopulation % 2 == 1) descendants.push_back(population[sizePopulation - 1]);
	return descendants;
}

Population EGA::PartialMappingCrossover() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> dist(0, numGenes - 1);

	for (int it = 0; it < sizePopulation / 2 * 2; it++) {
		Individ descendant(numGenes);
		Individ parentFir = population[it];
		Individ parentSec;
		if (it % 2 != 0) parentSec = population[it - 1];
		else parentSec = population[it + 1];

		int startPos = dist(gen);
		int endPos = dist(gen);
		if (endPos < startPos) {
			std::swap(startPos, endPos);
		}

		std::unordered_map<int, int> mapping;
		for (int i = startPos; i <= endPos; ++i) {
			descendant.encoding[i] = parentFir.encoding[i];
			mapping[parentFir.encoding[i]] = parentSec.encoding[i];
		}

		for (int i = 0; i < descendant.encoding.size(); ++i) {
			if (i >= startPos && i <= endPos) {
				continue;
			}
			int gene = parentSec.encoding[i];
			while (mapping.find(gene) != mapping.end()) {
				gene = mapping[gene];
			}
			descendant.encoding[i] = gene;
		}
		descendant.CalculateAdaptability(distances);
		descendants.push_back(descendant);
	}
	if (sizePopulation % 2 == 1) descendants.push_back(population[sizePopulation - 1]);
	return descendants;
}

Population EGA::PointMutation() {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(0.0, 1.0);
	std::uniform_int_distribution<int> gene(0, numGenes - 1);

	for (int it = 0; it < descendants.size(); it++) {
		Individ descendant = descendants[it];
		bool IsMutated = false;
		for (int i = 0; i < numGenes; i++) {
			if (dis(gen) < descendant.probabilityGeneMutation) {
				int geneToSwap;
				do {
					geneToSwap = gene(gen);
				} while (geneToSwap == i);
				descendant.SwapGenes(i, geneToSwap);
				IsMutated = true;
				break;
			}
		}
		if (IsMutated) {
			descendant.CalculateAdaptability(distances);
			mutants.push_back(descendant);
		}
	}
	return mutants;
}

Population EGA::StrategyParentsDescendantsMutants() {
	reproductionSet.insert(reproductionSet.end(), population.begin(), population.end());
	population.assign(sizePopulation, 0);
	reproductionSet.insert(reproductionSet.end(), descendants.begin(), descendants.end());
	descendants.clear();
	reproductionSet.insert(reproductionSet.end(), mutants.begin(), mutants.end());
	mutants.clear();
	return reproductionSet;
}

Population EGA::StrategyDescendantsMutants() {
	population.assign(sizePopulation, 0);
	reproductionSet.insert(reproductionSet.end(), descendants.begin(), descendants.end());
	descendants.clear();
	reproductionSet.insert(reproductionSet.end(), mutants.begin(), mutants.end());
	mutants.clear();
	return reproductionSet;
}

Population EGA::ProportionalSelection() {
	for (int i = 0; i < sizePopulation; i++) {
		std::vector<float> adaptabilityVector;
		std::transform(reproductionSet.begin(), reproductionSet.end(), std::back_inserter(adaptabilityVector),
			[](const Individ& individ) { return individ.adaptability; });

		Roulette rouletteIndivids = Roulette(adaptabilityVector);
		int indexIndivid = rouletteIndivids.ChooseObjectByRoulette();

		population[i] = reproductionSet[indexIndivid];
		reproductionSet.erase(reproductionSet.begin() + indexIndivid);
	}

	Individ newBestIndidvid = FindBestSolution();
	if (newBestIndidvid.adaptability > BestIndivid.adaptability) {
		BestIndivid = newBestIndidvid;
		NumIterationsWithoutImprovements = 0;
	}
	else NumIterationsWithoutImprovements++;

	reproductionSet.clear();
	return population;
}

Population EGA::RankSelection() {
	std::sort(reproductionSet.begin(), reproductionSet.end(), [](const Individ& a, const Individ& b) {
		return a.adaptability > b.adaptability;
		});

	std::vector<float> adaptabilities;
	float v = reproductionSet.size();
	for (int rank = 0; rank < v; rank++) {
		float value = 2 * (v - rank);
		adaptabilities.push_back(value);
	}

	for (int it = 0; it < sizePopulation; it++) {
		Roulette rouletteIndivids = Roulette(adaptabilities);
		int indexIndivid = rouletteIndivids.ChooseObjectByRoulette();

		population[it] = reproductionSet[indexIndivid];
		reproductionSet.erase(reproductionSet.begin() + indexIndivid);
		adaptabilities.erase(adaptabilities.begin() + indexIndivid);
	}

	Individ newBestIndidvid = FindBestSolution();
	if (newBestIndidvid.adaptability > BestIndivid.adaptability) {
		BestIndivid = newBestIndidvid;
		NumIterationsWithoutImprovements = 0;
	}
	else NumIterationsWithoutImprovements++;

	reproductionSet.clear();
	return population;
}

Population EGA::BetaTournamentSelection(int dimensionTournament) {
	for (int it = 0; it < sizePopulation; it++) {
		if (dimensionTournament > reproductionSet.size()) dimensionTournament = reproductionSet.size();

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<int> dist(0, reproductionSet.size() - 1);

		std::map<int, float> betaTournament;
		for (int i = 0; i < dimensionTournament; i++) {
			int indexIndivid = dist(gen);
			while (betaTournament.find(indexIndivid) != betaTournament.end()) indexIndivid = dist(gen);
			betaTournament.insert(std::make_pair(indexIndivid, reproductionSet[indexIndivid].adaptability));
		}

		auto iterator = std::max_element(betaTournament.begin(), betaTournament.end(),
			[](const std::pair<int, float>& a, const std::pair<int, float>& b) {
				return a.second < b.second;
			});

		if (iterator != betaTournament.end()) {
			Individ individToAdd = reproductionSet[iterator->first];
			population[it] = individToAdd;
			reproductionSet.erase(reproductionSet.begin() + iterator->first);
		}
	}

	Individ newBestIndidvid = FindBestSolution();
	if (newBestIndidvid.adaptability > BestIndivid.adaptability) {
		BestIndivid = newBestIndidvid;
		NumIterationsWithoutImprovements = 0;
	}
	else NumIterationsWithoutImprovements++;

	reproductionSet.clear();
	return population;
}

bool EGA::StopConditionNumIterations() {
	if (++numIterations == maxNumIterations) return true;
	return false;
}

Individ EGA::FindBestSolution() {
	return *std::max_element(population.begin(), population.end(), [](const Individ& a, const Individ& b) {
		return a.adaptability < b.adaptability;
		});
}

bool EGA::StopConditionNoImprovements() {
	if (NumIterationsWithoutImprovements == maxNumIterationsWithoutImprovements) return true;
	return false;
}

void EGA::PrintResultsPopulation() {
	std::cout << "Особи популяции:" << std::endl;
	for (int i = 0; i < population.size(); i++) {
		Individ curIndivid = population[i];
		std::cout << "Кодировка особи - ";
		for (int j = 0; j < curIndivid.encoding.size(); j++) {
			std::cout << curIndivid.encoding[j] << " ";
		}
		std::cout << "   Приспособленность особи - " << curIndivid.adaptability << std::endl;
	}
	std::cout << std::endl << "Лучшая особь - ";
	for (int j = 0; j < BestIndivid.encoding.size(); j++) {
		std::cout << BestIndivid.encoding[j] << " ";
	}
	std::cout << "   Приспособленность особи - " << BestIndivid.adaptability << std::endl << std::endl;
}
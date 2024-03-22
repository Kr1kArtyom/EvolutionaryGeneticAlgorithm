#pragma once
#include <vector> 
#include <random>
#include <algorithm>

typedef std::vector<std::vector<float>> MatrixDistances;

struct Roulette {
    std::vector<float> roulette;
    float sizeRoulette;

    Roulette(std::vector<float> adaptabilities);
    int ChooseObjectByRoulette();
};

struct Individ {
    std::vector<int> encoding;
    float adaptability;
    float probabilityGeneMutation;

    Individ(int numCities = 0, float _probabilityGeneMutation = 0.01);
    float CalculateAdaptability(const MatrixDistances distances);
    void SwapGenes(int indFirst, int indSecond);
};

typedef std::vector<Individ> Population;

class EGA
{
private:
    int numGenes, sizePopulation;
    MatrixDistances distances;
    Population population;
    Population descendants;
    Population mutants;
    Population reproductionSet;

    int maxNumIterations, numIterations;
    Individ BestIndivid;
    int maxNumIterationsWithoutImprovements, NumIterationsWithoutImprovements;
public:
    EGA(int _sizePopulation, MatrixDistances _distances, int _maxNumIterations = 10, int _maxNumIterationsWithoutImprovements = 5);
    void swapIndivids(Individ& indFirst, Individ& indSecond);
    Population generateRandomPopulations();
    Population generateInitialPopulation();
    Population generateGreedyPopulations();

    Population RandomSelection();
    Population PositiveAssortativeMating();
    Population NegativeAssortativeMating();
    Population InbreedingSelection();
    Population OutbreedingSelection();

    Population OrderCrossover();
    Population PartialMappingCrossover();

    Population PointMutation();

    Population StrategyParentsDescendantsMutants();
    Population StrategyDescendantsMutants();

    Population ProportionalSelection();
    Population RankSelection();
    Population BetaTournamentSelection(int dimensionTournament = 3);

    bool StopConditionNumIterations();
    bool StopConditionNoImprovements();

    Individ FindBestSolution();

    void PrintResultsPopulation();
};
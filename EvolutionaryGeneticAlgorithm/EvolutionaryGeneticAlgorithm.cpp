#include <iostream>
#include "EGA.h"

using namespace std;

void Menu(vector<vector<float>>& distances) {
	int sizePopulation;
	cout << "Введите размер популяции: ";
	cin >> sizePopulation;

	EGA ega(sizePopulation, distances);

	int initialPopulation;
	cout << endl << "Выберите оператор формирования начальной популяции" << endl
		<< "1) Случайно" << endl
		<< "1) Жадный алгоритм - Метод ближайшего города" << endl
		<< "2) Жадный алгоритм - Метод ближайшего соседа" << endl;
	cin >> initialPopulation;

	int parentPairs;
	cout << endl << "Выберите оператор выбора родительской пары" << endl
		<< "1) Случайно" << endl
		<< "2) Инбридинг" << endl
		<< "3) Аутбридинг" << endl
		<< "4) Положительное ассоциативное скрещивание" << endl
		<< "5) Отрицательное ассоциативное скрещивание" << endl;
	cin >> parentPairs;

	int crossover;
	cout << endl << "Веберите оператор кроссовера" << endl
		<< "1) Порядковый" << endl
		<< "2) Частичное отображение" << endl;
	cin >> crossover;

	int mutation;
	cout << endl << "Выберите опертор мутации" << endl
		<< "1) Точечная" << endl;
	cin >> mutation;

	int strategyFormationNextPop;
	cout << endl << "Выберите стратегию формирования следующей популяции" << endl
		<< "1) Родители + потомки + мутанты" << endl
		<< "2) Потомки + мутанты" << endl;
	cin >> strategyFormationNextPop;

	int selection;
	cout << endl << "Выберите оператор селекции" << endl
		<< "1) Пропорциональная" << endl
		<< "2) Ранговая" << endl
		<< "3) beta - турнир" << endl;
	cin >> selection;

	switch (initialPopulation) {
	case 1:
		ega.generateRandomPopulations();
		break;
	case 2:
		ega.generateGreedyPopulations();
		break;
	case 3:
		ega.generateInitialPopulation();
	}

	do {
		switch (parentPairs) {
		case 1:
			ega.RandomSelection();
			break;
		case 2:
			ega.InbreedingSelection();
			break;
		case 3:
			ega.OutbreedingSelection();
			break;
		case 4:
			ega.PositiveAssortativeMating();
			break;
		case 5:
			ega.NegativeAssortativeMating();
			break;
		}

		switch (crossover) {
		case 1:
			ega.OrderCrossover();
			break;
		case 2:
			ega.PartialMappingCrossover();
			break;
		}

		switch (mutation) {
		case 1:
			ega.PointMutation();
			break;
		}

		switch (strategyFormationNextPop) {
		case 1:
			ega.StrategyParentsDescendantsMutants();
			break;
		case 2:
			ega.StrategyDescendantsMutants();
			break;
		}

		switch (selection) {
		case 1:
			ega.ProportionalSelection();
			break;
		case 2:
			ega.RankSelection();
			break;
		case 3:
			ega.BetaTournamentSelection();
			break;
		}

		ega.PrintResultsPopulation();
	} while (!ega.StopConditionNoImprovements() && !ega.StopConditionNumIterations());
}

int main()
{
	setlocale(LC_ALL, "rus");

	vector<vector<float>> distances{
		{0.00, 12.53, 5.10, 11.70, 3.00, 13.04, 4.47, 1.00, 8.49, 7.62, 10.00, 13.15, 13.00, 13.60, 13.15},
		{ 12.53, 0.00, 10.05, 10.00, 11.40, 7.28, 13.45, 13.04, 13.00, 8.06, 5.39, 8.25, 6.32, 2.83, 4.47 },
		{ 5.10, 10.05, 0.00, 13.45, 2.24, 13.42, 9.06, 6.08, 12.08, 2.83, 5.83, 13.89, 13.00, 12.04, 12.37 },
		{ 11.70, 10.00, 13.45, 0.00, 13.04, 3.61, 9.00, 11.40, 5.39, 13.60, 13.00, 2.83, 4.47, 8.25, 6.32 },
		{ 3.00, 11.40, 2.24, 13.04, 0.00, 13.60, 7.28, 4.00, 10.82, 5.00, 7.81, 13.93, 13.34, 13.04, 13.04 },
		{ 13.04, 7.28, 13.42, 3.61, 13.60, 0.00, 11.40, 13.00, 8.60, 12.81, 11.40, 1.00, 1.00, 5.00, 3.00 },
		{ 4.47, 13.45, 9.06, 9.00, 7.28, 11.40, 0.00, 3.61, 4.47, 11.05, 12.65, 11.18, 11.70, 13.60, 12.53 },
		{ 1.00, 13.04, 6.08, 11.40, 4.00, 13.00, 3.61, 0.00, 7.81, 8.54, 10.82, 13.04, 13.04, 13.93, 13.34 },
		{ 8.49, 13.00, 12.08, 5.39, 10.82, 8.60, 4.47, 7.81, 0.00, 13.34, 14.00, 8.06, 9.22, 12.21, 10.63 },
		{ 7.62, 8.06, 2.83, 13.60, 5.00, 12.81, 11.05, 8.54, 13.34, 0.00, 3.16, 13.45, 12.21, 10.44, 11.18 },
		{ 10.00, 5.39, 5.83, 13.00, 7.81, 11.40, 12.65, 10.82, 14.00, 3.16, 0.00, 12.21, 10.63, 8.06, 9.22 },
		{ 13.15, 8.25, 13.89, 2.83, 13.93, 1.00, 11.18, 13.04, 8.06, 13.45, 12.21, 0.00, 2.00, 6.00, 4.00 },
		{ 13.00, 6.32, 13.00, 4.47, 13.34, 1.00, 11.70, 13.04, 9.22, 12.21, 10.63, 2.00, 0.00, 4.00, 2.00 },
		{ 13.60, 2.83, 12.04, 8.25, 13.04, 5.00, 13.60, 13.93, 12.21, 10.44, 8.06, 6.00, 4.00, 0.00, 2.00 },
		{ 13.15, 4.47, 12.37, 6.32, 13.04, 3.00, 12.53, 13.34, 10.63, 11.18, 9.22, 4.00, 2.00, 2.00, 0.00 }
	};

	Menu(distances);
}
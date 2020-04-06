/**
Calculate stabilometry index from jiang_CoPpositon_pos.sto
based on this article:
Kaoru Imaoka, Hitoshi Murase, Miho Fukuhara: "Collection of Data for Healthy Subjects in Stabilometry," Equilibrium Res Suppl. Vol 12, pp. 1-84, 1997.
*/
#include <OpenSim/OpenSim.h>
#include <vector>
#include <iostream>
#include <fstream>

using namespace OpenSim;
using namespace SimTK;
using namespace std;

double getLengthIn1D(const vector<double>& data)
{
	if (data.size() <= 1) return 0.0;

	vector<double> diffAbsVec(data.size() - 1, 0.0);
	for (int i = 1; i < data.size(); i++) {
		diffAbsVec.push_back(abs(data[i] - data[i - 1L]));
	}
	// Sum up in ascending order to avoid information loss error
	sort(diffAbsVec.begin(), diffAbsVec.end());
	double sum = 0.0;
	for (double value : diffAbsVec) {
		sum += value;
	}
	return sum;
}

double getLengthIn2D(const vector<double>& xdata, const vector<double>& ydata)
{
	if (xdata.size() <= 1 || ydata.size() <= 1 || xdata.size() != ydata.size()) return 0.0;
	vector<double> deltaLengthVec(xdata.size() - 1, 0.0);
	for (int i = 1; i < xdata.size(); i++) {
		deltaLengthVec.push_back(sqrt((xdata[i] - xdata[i - 1L]) * (xdata[i] - xdata[i - 1L]) + (ydata[i] - ydata[i - 1L]) * (ydata[i] - ydata[i - 1L])));
	}
	// Sum up in ascending order to avoid information loss error
	sort(deltaLengthVec.begin(), deltaLengthVec.end());
	double sum = 0.0;
	for (double value : deltaLengthVec) {
		sum += value;
	}
	return sum;
}

double getAverage(const vector<double>& data)
{
	return accumulate(data.begin(), data.end(), 0.0) / data.size();
}

double getStandardDeviation(const vector<double>& data)
{
	const double mean = getAverage(data);
	vector<double> devVec(data.size());
	for (double x : data) {
		devVec.push_back((x - mean) * (x - mean));
	}
	// Sum up in ascending order to avoid information loss error
	sort(devVec.begin(), devVec.end());
	double sum = 0.0;
	for (double value : devVec) {
		sum += value;
	}
	return sqrt(sum / data.size());
}

double getAverageOneWayVelocity(bool isPositiveDirection, const vector<double>& data, double DT)
{
    if (data.size() <= 1) return 0.0;

    vector<double> deltaVec(0);
    if (isPositiveDirection) {
        for (int i = 1; i < data.size(); ++i) {
            const double delta = data[i] - data[i-1];
            if (delta > 0.0) deltaVec.push_back(delta);
        }
    } else {
        for (int i = 1; i < data.size(); ++i) {
            const double delta = data[i] - data[i-1];
            if (delta < 0.0) deltaVec.push_back(-delta);
        }
    }
    // Sum up in ascending order to avoid information loss error
    sort(deltaVec.begin(), deltaVec.end());
    return getAverage(deltaVec)/DT;
}

double getRectangleArea(const vector<double>& xdata, const vector<double>& ydata)
{
	const double xMin = *min_element(xdata.begin(), xdata.end());
	const double xMax = *max_element(xdata.begin(), xdata.end());
	const double yMin = *min_element(ydata.begin(), ydata.end());
	const double yMax = *max_element(ydata.begin(), ydata.end());
	return (xMax - xMin) * (yMax - yMin);
}

double getSurroundingArea(vector<double> xdata, vector<double> ydata)
{
	const int size = min(xdata.size(), ydata.size());
	const double xMean = getAverage(xdata);
	const double yMean = getAverage(ydata);
	// Translation around average point
	for (int i = 0; i < size; ++i) {
		xdata[i] -= xMean;
		ydata[i] -= yMean;
	}
	// Polar coordinate transformation
	vector<double> rData(size);
	vector<double> thetaData(size); // in radian
	for (int i = 0; i < size; ++i) {
		const double x = xdata[i];
		const double y = ydata[i];
		rData[i] = sqrt(x * x + y * y);
		thetaData[i] = atan2(y, x); // return in [-pi,pi]
	} // TODO: sort (r,theta) by theta
	// Extract max radius in the angle(deltaTheta)
	const int deltaTheta = 3; // in degree
	const int divisionNumber = 360 / deltaTheta;
	vector<double> rMaxVec(divisionNumber, 0.0);
	for (int i = 0; i < rMaxVec.size(); ++i) {
		const double thetaMin = (i * deltaTheta - 180) * SimTK_DEGREE_TO_RADIAN;
		const double thetaMax = ((i + 1) * deltaTheta - 180) * SimTK_DEGREE_TO_RADIAN;
		for (int j = 0; j < size; ++j) {
			const double theta = thetaData[j];
			if (theta > thetaMin&& theta <= thetaMax && rData[j] > rMaxVec[i]) {
				rMaxVec[i] = rData[j];
			}
		}
	}
	vector<double> subAreaVec;
	int subDivisionCount = 0;
	int thetaStartIndex = 0;
	while (rMaxVec[thetaStartIndex % divisionNumber] <= 0.0) ++thetaStartIndex;
	while (subDivisionCount < divisionNumber) {
		int thetaEndIndex = thetaStartIndex + 1;
		while (rMaxVec[thetaEndIndex % divisionNumber] <= 0.0) ++thetaEndIndex;
		subAreaVec.push_back(0.5 * rMaxVec[thetaStartIndex % divisionNumber] * rMaxVec[thetaEndIndex % divisionNumber] * sin((thetaEndIndex - thetaStartIndex) * deltaTheta * SimTK_DEGREE_TO_RADIAN));
		subDivisionCount += thetaEndIndex - thetaStartIndex;
		thetaStartIndex = thetaEndIndex;
	}
	if (subDivisionCount != divisionNumber) cout << "Divison error." << endl;
	// Sum up in ascending order to avoid information loss error
	sort(subAreaVec.begin(), subAreaVec.end());
	double sum = 0.0;
	for (double value : subAreaVec) {
		sum += value;
	}
	return sum;
}

int main()
{
	Storage sto("jiang_CoPpositon_pos.sto");
    
    const double DT = 1.0/20;
    //sto.resample(DT, 3);
    sto.resampleLinear(DT); // set the sampling rate

    sto.crop(2.0, 4.0);
    
	const int timeSize = sto.getSize();
	cout << "timeSize: " << timeSize << endl;
	const double timeRange = sto.getLastTime() - sto.getFirstTime();
	cout << "timeRange: " << timeRange << " s" << endl;

	// Get CoP-x and CoP-z data
	Array<double> CoPxArray, CoPzArray;
	sto.getDataColumn("cop-x", CoPxArray);
	sto.getDataColumn("cop-z", CoPzArray);
	/*
	for (int i = 0; i < CoPxArray.getSize();i++) {
		cout << i << "\tx: " << CoPxArray[i] << "\ty: " << CoPzArray[i] << endl;
	}
	*/
	// Convert OpenSim::Array to std::vector
	std::vector<double> CoPxVector(CoPxArray.getSize());
	std::vector<double> CoPzVector(CoPzArray.getSize());
	for (int i = 0; i < timeSize; ++i) {
		CoPxVector[i] = CoPxArray[i];
		CoPzVector[i] = CoPzArray[i];
	}

	const double totalLength = getLengthIn2D(CoPzVector, CoPxVector);
	cout << "Total length: " << totalLength << " m" << endl;
	const double totalLengthPerUnitTime = totalLength / timeRange;
	cout << "Total length per unit time: " << totalLengthPerUnitTime << " m/s" << endl;

	const double surroundingArea = getSurroundingArea(CoPzVector, CoPxVector);
	const double totalLengthPerUnitArea = totalLength / surroundingArea;
	cout << "Total length per unit area: " << totalLengthPerUnitArea << " m^(-1)" << endl;

	const double APlength = getLengthIn1D(CoPxVector);
	const double MLlength = getLengthIn1D(CoPzVector);
	cout << "AP length: " << APlength << " m" << endl;
	cout << "ML length: " << MLlength << " m" << endl;
	const double APlengthPerUnitTime = APlength / timeRange;
	const double MLlengthPerUnitTIme = MLlength / timeRange;
	cout << "AP length per unit time: " << APlengthPerUnitTime << " m/s" << endl;
	cout << "ML length per unit time: " << MLlengthPerUnitTIme << " m/s" << endl;

	const double rectangleArea = getRectangleArea(CoPzVector, CoPxVector);
	cout << "Rectangle area: " << rectangleArea << " m^2 (" << rectangleArea * 1.0e+4 << " cm^2)" << endl;
	cout << "Surrounding area: " << surroundingArea << " m^2 (" << surroundingArea * 1.0e+4 << " cm^2)" << endl;

	const double APSD = getStandardDeviation(CoPxVector);
	const double MLSD = getStandardDeviation(CoPzVector);
	cout << "AP SD: " << APSD << " m" << endl;
	cout << "ML SD: " << MLSD << " m" << endl;
    
    const double AnteriorVel = getAverageOneWayVelocity(true, CoPxVector, DT);
    const double PosteriorVel = getAverageOneWayVelocity(false, CoPxVector, DT);
    cout << "Average Anterior Velocity: " << AnteriorVel << " m/s" << endl;
    cout << "Average Posterior Velocity: " << PosteriorVel << " m/s" << endl;

	ofstream outfile("results.csv");
	outfile << "total length/m,total length per unit time/(m/s),total length per unit area//m,AP length/m,ML length/m,AP length per unit time/(m/s),ML length per unit time/(m/s),rectangle area/m2,surrounding area/m2,AP SD/m,ML SD/m,average anterior velocity/(m/s),average posterior velocity/(m/s)" << endl;
	outfile << totalLength << "," << totalLengthPerUnitTime << "," << totalLengthPerUnitArea << "," << APlength << "," << MLlength << "," << APlengthPerUnitTime << "," << MLlengthPerUnitTIme << "," << rectangleArea << "," << surroundingArea << "," << APSD << "," << MLSD << "," << AnteriorVel << "," << PosteriorVel << endl;
	outfile.close();

	return 0;
}

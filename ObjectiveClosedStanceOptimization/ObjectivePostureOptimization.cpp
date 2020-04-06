/* --------------------------------------------------------------------------*
* This program was modified from OpenSim: SimpleOptimizationExample.cpp      *
*                                                                            *
* Original copyright (c) 2005-2017 Stanford University and the Authors       *
* Original author(s): Ayman Habib                                            *
*                                                                            *
* Licensed under the Apache License, Version 2.0 (the "License"); you may    *
* not use this file except in compliance with the License. You may obtain a  *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
*                                                                            *
* Unless required by applicable law or agreed to in writing, software        *
* distributed under the License is distributed on an "AS IS" BASIS,          *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
* See the License for the specific language governing permissions and        *
* limitations under the License.                                             *
* -------------------------------------------------------------------------- */

//==============================================================================
//==============================================================================
#include <OpenSim/OpenSim.h>
#include <ctime>  // clock(), clock_t, CLOCKS_PER_SEC
#include <iostream>
#include <cstdio>
#include <mutex>
#include <fstream>

using namespace OpenSim;
using namespace SimTK;
using namespace std;

// step count for troubleshooting
int stepCount = 0;

double bestSoFar = Infinity;
State bestState;

mutex mtx_objectiveFunc;

class ObjectivePostureOptimizationSystem : public OptimizerSystem {
public:
	/* Constructor class. Parameters passed are accessed in the objectiveFunc() class. */
	ObjectivePostureOptimizationSystem(int numParameters, State& s, Model& aModel) :
		OptimizerSystem(numParameters),
		si(s),
		osimModel(aModel) {}

	int objectiveFunc(const Vector& newControls, bool new_coefficients, Real& f) const override {
		// make a copy of out initial states
		State s = si;

		// Set values of coordinates
		const CoordinateSet& coords = osimModel.getCoordinateSet();
		coords.get(7).setValue(s, newControls[0]);
		for (int i = 0; i < 6; i++) {
			coords.get(i + 9).setValue(s, newControls[i + 1]);
			coords.get(i + 15).setValue(s, newControls[i + 1]);
		}
		coords.get(21).setValue(s, newControls[7]);

		const Set<Muscle>& muscleSet = osimModel.getMuscles();
		// Ignore muscle force
		for (int i = 0; i < muscleSet.getSize(); i++) {
			muscleSet.get(i).setIgnoreActivationDynamics(s, true);
			//muscleSet.get(i).setAppliesForce(s, false);
		}
		// Make sure the muscles states are in equilibrium
		osimModel.equilibrateMuscles(s);

		// Realize the state to Dynamics stage
		if (s.getSystemStage() < SimTK::Stage::Dynamics) {
			osimModel.realizeDynamics(s);
		}

		// Calculate Center of Mass
		Vec3 com(0.0);
		double total_mass = 0.0;
		const BodySet& bodies = osimModel.getBodySet();
		for (int i = 0; i < bodies.getSize(); i++) {
			const OpenSim::Body& body = bodies[i];
			if ("platform" == body.getName() || "ground" == body.getName()) {
				continue;
			}
			total_mass += body.getMass();
			com += body.getMass() * body.findStationLocationInGround(s, body.getMassCenter());
		}
		com /= total_mass;

		double J_static = 0.0;
		double bos_center_x = 0.0;
		double cop_x = 0.0;
		Array<double> rcf, lcf; // right / left calcn force
		rcf = osimModel.getForceSet().get("foot_floor_r").getRecordValues(s);
		lcf = osimModel.getForceSet().get("foot_floor_l").getRecordValues(s);
		if (rcf[1] < 0.0 && lcf[1] < 0.0) { // reaction force < 0 -> the person is grounding
		  // Calculate J_static
			const Vec3 heel_r_position_in_ground = osimModel.getContactGeometrySet().get("foot_r").getFrame().findStationLocationInGround(s, Vec3(0.0));
			bos_center_x = heel_r_position_in_ground[0] + 0.10127488;
			cop_x = (rcf[5] + lcf[5]) / (rcf[1] + lcf[1]); // Convert torques and forces to CoP (M = r x F where r = {CoPx, 0, CoPz})
			const double com_x = com[0];
			J_static = (cop_x - bos_center_x) * (cop_x - bos_center_x) + (com_x - bos_center_x) * (com_x - bos_center_x);
		}
		else { // if the person is floating
			// add penalty
			J_static = fabs(com[1]); // height of Center of Mass
		}

		// Calculate J_omega
		double J_omega = 0.0;
		Vector torques;
		osimModel.getMatterSubsystem().multiplyBySystemJacobianTranspose(s, osimModel.getMultibodySystem().getRigidBodyForces(s, SimTK::Stage::Dynamics), torques);
		/**
		torques[0]: platform_rz, torques[1]: platform_tx, torques[2]: platform_tz,
		torques[3]: pelvis_tilt, torques[4]: pelvis_list,
		torques[5]: pelvis_rotation, torques[6]: pelvis_tx, torques[7]: pelvis_ty, torques[8]: pelvis_tz,
		torques[9]: hip_flexion_r, torques[10]:hip_adduction_r, torques[11]:hip_rotation_r,
		torques[12]:hip_flexion_l, torques[13]:hip_adduction_l, torques[14]:hip_rotation_l,
		torques[15]:lumbar_extension, torques[16]:lumbar_bending, torques[17]:lumbar_rotation,
		torques[18]:knee_angle_r,
		torques[19]:knee_angle_l,
		torques[20]:ankle_angle_r,
		torques[21]:ankle_angle_l,
		torques[22]:subtalar_angle_r,
		torques[23]:subtalar_angle_l
		*/
		for (int i = 9; i < 24; i++) {
			if (0.0 == torques[i]) { // if the posture is invalid
				J_omega += 1.0e+3;
			}
			else {
				J_omega += torques[i] * torques[i];
			}
		}

		// Calculate J_close
		double J_close = 0.0;
		const MarkerSet& markerSet = osimModel.getMarkerSet();
		/*
		const double distance_metatarsal = markerSet.get("metatarsal_r").calcDistanceBetween(s, markerSet.get("metatarsal_l"));
		const double distance_heel = markerSet.get("heel_r").calcDistanceBetween(s, markerSet.get("heel_l"));
		*/
		const double distance_metatarsal = markerSet.get("metatarsal_r").getLocationInGround(s)[2] - markerSet.get("metatarsal_l").getLocationInGround(s)[2];
		const double distance_heel = markerSet.get("heel_r").getLocationInGround(s)[2] - markerSet.get("heel_l").getLocationInGround(s)[2];
		J_close = distance_metatarsal + distance_heel;
		if (J_close < 0.0) J_close = -10.0 * J_close; // if the foots are crossed

		// Calculate objective function
		const double w_static = 1.0e+4;
		const double w_omega = 1.0e-1;
		const double w_close = 1.0e+2;
		f = w_static * J_static + w_omega * J_omega + w_close * J_close;
		const double normal_reaction_force = rcf[7] + lcf[7];
		f += normal_reaction_force; // Add normal reaction force to decrease intrusion to floor

		{
			lock_guard<mutex> lock(mtx_objectiveFunc);
			stepCount++;
			if (f < bestSoFar) {
				bestSoFar = f;
				bestState = s;
				cout << "\nobjective evaluation #: " << stepCount << endl;
				for (int i = 0; i < newControls.size(); i++) {
					cout << "controls[" << i << "] " << newControls[i] << endl; // print controls in degree
				}
				for (int i = 0; i < 9; i++) {
					cout << i << "\t" << osimModel.getCoordinateSet().get(i).getName() << "\t" << osimModel.getCoordinateSet().get(i).getLocked(s) << "\t" << osimModel.getCoordinateSet().get(i).getValue(s) << endl;
				}
				for (int i = 9; i < osimModel.getNumCoordinates(); i++) {
					cout << i << "\t" << osimModel.getCoordinateSet().get(i).getName() << "\t" << osimModel.getCoordinateSet().get(i).getLocked(s) << "\t" << osimModel.getCoordinateSet().get(i).getValue(s) * SimTK_RADIAN_TO_DEGREE << endl;
				}
				for (int i = 0; i < rcf.size(); i++) {
					cout << "rcf[" << i << "]\t" << osimModel.getForceSet().get("foot_floor_r").getRecordLabels()[i] << "\t" << rcf[i] << endl;
				}
				for (int i = 0; i < lcf.size(); i++) {
					cout << "lcf[" << i << "]\t" << osimModel.getForceSet().get("foot_floor_l").getRecordLabels()[i] << "\t" << lcf[i] << endl;
				}
				for (int i = 0; i < torques.size(); i++) {
					cout << "torques[" << i << "] " << torques[i] << endl;
				}
				cout << "w*Jstatic: " << w_static * J_static << "\tJstatic: " << J_static << endl;
				cout << "w*Jomega: " << w_omega * J_omega << "\tJomega: " << J_omega << endl;
				cout << "w*Jclose: " << w_close * J_close << "\tJclose: " << J_close << endl;
				cout << "bos: " << bos_center_x << endl;
				cout << "cop_x: " << cop_x << endl;
				cout << "com: " << com << endl;
				cout << "f: " << f << endl;
			}
		}

		return(0);
	}

private:
	State& si;
	Model& osimModel;
};

int main()
{
	try {
		std::clock_t startTime = std::clock();

		// Create a new OpenSim model
		Model osimModel("new_model_03.osim");

		// Make sure of ElasticFoundationForce contact between foot and floor so as to calculate CoP
		if (osimModel.getForceSet().get("foot_floor_r").getConcreteClassName() != "ElasticFoundationForce" || osimModel.getForceSet().get("foot_floor_l").getConcreteClassName() != "ElasticFoundationForce") {
			cout << "Contacts of foot and floor are not ElasticFoundationForce!" << endl;
			return 1;
		}

		// Load markers from file
		MarkerSet markers("markers.xml");
		osimModel.updateMarkerSet(markers);

		// Initialize the system and get the state representing the state system
		State& si = osimModel.initSystem();

		const CoordinateSet& coords = osimModel.getCoordinateSet();

		// Initialize the optimizer system we've defined.
		const int num_controls = 8;
		ObjectivePostureOptimizationSystem sys(num_controls, si, osimModel);
		Real f = NaN;

		/* Define initial values and bounds for the controls to optimize */
		Vector controls(num_controls, 0.0); // 0 radian for default value
		Vector lower_bounds(num_controls, 0.0);
		Vector upper_bounds(num_controls, 0.0);

		// Define initial values
		controls[0] = 0.94; // pelvis_ty
		controls[1] = -5.0; // hip_flexion
		controls[2] = 3.0; // hip_adduction
		controls[3] = 2.0; // hip_rotation
		controls[4] = 1.0; // knee_angle
		controls[5] = 2.0; // ankle_angle
		controls[6] = 4.0; // subtalar_angle
		controls[7] = -10.0; // lumbar_extension
		// Convert degrees to radians
		Units degrees(Units::UnitType::Degrees);
		for (int i = 1; i < num_controls; i++) {
			controls[i] = degrees.convertTo(Units::UnitType::Radians, controls[i]);
		}

		// Define bounds
		lower_bounds[0] = 0.91; // coords[7].getRangeMin(); // Decrease intrusion to floor
		upper_bounds[0] = 0.97; // coords[7].getRangeMax();
		const double five_degree_in_radian = degrees.convertTo(Units::UnitType::Radians, 5.0);
		for (int i = 1; i < num_controls; i++) {
			lower_bounds[i] = controls[i] - five_degree_in_radian;
			upper_bounds[i] = controls[i] + five_degree_in_radian;
		}
		sys.setParameterLimits(lower_bounds, upper_bounds);

		for (int i = 0; i < num_controls; i++) {
			if (0 == i) {
				cout << i << "\t" << lower_bounds[i] << "\t" << upper_bounds[i] << endl;
			}
			else {
				cout << i << "\t" << lower_bounds[i] * SimTK_RADIAN_TO_DEGREE << "\t" << upper_bounds[i] * SimTK_RADIAN_TO_DEGREE << endl;
			}
		}

		// Create an optimizer. Pass in our OptimizerSystem
		// and the name of the optimization algorithm.
		Optimizer opt(sys, SimTK::CMAES);

		// Specify settings for the optimizer
		opt.setConvergenceTolerance(0.001);
		opt.setMaxIterations(10000);
		opt.setDiagnosticsLevel(3); // Set the level of debugging info displayed.
		//opt.setAdvancedIntOption("seed", 4985); // Set the seed to generate identical results
		opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);
		opt.setAdvancedStrOption("parallel", "multithreading");
		int nthreads = ParallelExecutor::getNumProcessors();
		cout << "Number of threads: " << nthreads << endl;
		//opt.setAdvancedIntOption("nthreads", 10);

		cout << "Ready to optimize." << endl;

		// Optimize it!
		f = opt.optimize(controls);
		cout << "Best f = " << f << endl;
		cout << "Elapsed time = " << (std::clock() - startTime) / CLOCKS_PER_SEC << "s" << endl;
		cout << "Objective posture optimization completed successfully." << endl;

		// Output the solution
		Manager manager(osimModel, bestState);
		manager.setSessionName(to_string(f));
		manager.integrate(0.0);
		Storage storage = manager.getStateStorage();
		storage.resampleLinear(1.0e-3);
		storage.print(to_string(f) + ".sto");
	}
	catch (const std::exception & ex)
	{
		std::cout << ex.what() << std::endl;
		return 1;
	}

	getchar();
	// End of main() routine.
	return 0;
}

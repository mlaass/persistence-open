/*
The main evaluation system for SDS persistence.
Moritz Laass, Martin Werner

*/
#define OVERWRITE_PROGRESS 1
// Default Headers for C++
#include <vector>
#include <cmath>
#include <set>
#include <numeric>
#include <limits>
#include <functional>
#include <sstream>
#include <chrono>

#include <fstream>
#include <iostream>
//#include <boost/range/irange.hpp>

// Boost Geometry
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

namespace bg = boost::geometry;
typedef bg::model::point<double, 2, bg::cs::cartesian> point;
typedef bg::model::linestring<point> linestring;
const auto as_str = [](const linestring &ls) {std::stringstream ss; ss << bg::wkt(ls); return ss.str(); };

#include "trajcomp-persistence.hpp"
#include <frechet/frechetrange.hpp>
#include "frechet_boostgeom.hpp"

// Diverse Trajectory Algorithm Implementations
#include <trajcomp/trajcomp.hpp>

// Configuration and aggregate result communication based on JSON
#include <jsonconf.hpp>

// Full result collection as CSV using emit lambda
typedef std::function<void(std::vector<std::string>)> emit_function_type;

/*
Section 1: Simplification Functions

This section implements relevant simplification functions for trajectories

- douglas_peucker(ls,eps):
       Douglas-Peucker with distance threshold eps
- persistence_beta(ls, beta):
       Persistence with Beta-Pruning, beta is the length of the bars to prune
- persistence_mrs(ls, beta, eps, rounds)
       Persistence simplificaiton with beta-pruning followed by <rounds> MRS rounds using
       a point distance threshold of eps. This threshold grows exponentially with rounds.
- persistence_sds(ls, beta, eps, rounds)
       Persistence simplificaiton with beta-pruning followed by <rounds> SDS rounds using
       a segment distance threshold of eps. This threshold grows exponentially with rounds.

*/

linestring douglas_peucker(const linestring &ls, double eps)
{
	auto traj = extract(ls);
	auto simplified = trajcomp::douglas_peucker(traj, eps);
	linestring ret;
	for (const auto &c : simplified)
		ret.push_back(point(c[0], c[1]));
	return ret;
}

linestring persistence_beta(const linestring &ls, double beta)
{
	return extract_linestring(persistence::persistence(persistence::traj_to_curve(extract(ls)), beta));
}

linestring persistence_mrs(const linestring &ls, double beta, double eps, int levels)
{
	auto traj = extract(ls);
	auto curve = persistence::traj_to_curve(traj);
	for (size_t i = 0; i < levels; i++)
	{
		auto p_result = persistence::persistence(curve, beta, -1);
		curve = prune_curve_dist(p_result.pruned, eps * pow(2, i));
	}
	linestring rls;
	for (auto c : curve)
		rls.push_back(ls[c.index]);
	return rls;
}

linestring persistence_sds(const linestring &ls, double beta, double eps, int levels)
{
	auto traj = extract(ls);
	auto curve = persistence::traj_to_curve(traj);
	for (size_t i = 0; i < levels; i++)
	{
		auto p_result = persistence::persistence(curve, beta, -1);
		//curve = prune_curve_dist(, pow(2,i));
		curve = prune_sds(p_result.pruned, eps);
		// use curve.index
		eps *= 2; // exponential growth as well - otherwise unfair
	}
	linestring rls;
	for (auto c : curve)
		rls.push_back(ls[c.index]);
	return rls;
}

JSONConfig cfg;

/*
Section 2: Algorithms for processing WKT files.

Algorithm 2.1: Scan a WKT file

scan (input, output, func):
   scan over input line by line, parse as WKT ignoring errors and call func for
   each line. Func signature is
       void func(ls, emit, first),
   where ls is the linestring, emit is a lambda writing vectors of strings to output and
   first is set in case of the first line is to be written (e.g., CSV headers)

*/

template <typename callback>
void scan(std::string input, std::string output, callback func)
{
	std::ifstream ifs(input);
	std::ofstream ofs(output);
	std::string line;
	if (ifs.fail())
		std::cout << "failed to find: " << input << std::endl;

	if (ofs.fail())
		std::cout << "failed to find: " << input << std::endl;

	ofs << std::setprecision(10);
	auto emit = [&ofs](const std::vector<std::string> &columns) {
		std::stringstream ss;
		std::copy(columns.begin(), columns.end(), std::ostream_iterator<std::string>(ss, ";"));
		std::string line = ss.str();
		line.pop_back();
		ofs << line << std::endl;
	};
	bool first = true;
	while (std::getline(ifs, line))
	{
		linestring ls;
		try
		{
			bg::read_wkt(line, ls);
		}
		catch (...)
		{
			//std::cout << "Not a geom: " << line << std::endl;
			continue;
		}
		// tune DP
		func(ls, emit, first);
		first = false;
	}
}

/*
Algorithm 2.2: Calibrate a simplification algorithm

calibrate (filename, N, s, interval, [tag])
   Determines good parameters for a single-parameter simplifier s reaching an average
   trajectory size of N on data taken from filename. The single parameter is initialized from
   an <interval>. As all our parameters are monotonous, we apply exponential growth to
   the second bound followed by nested interval search. Progress messages are decorated
   with the given tag to know which algorithm is running and how far the interval converged

   Data is written to "garbage.tmp" as scan does not yet allow to scan without output.

   Reacts to #define OVERWRITE_PROGRESS by either using carriage return and flush for
   interactive use  or simple newlines for log-like files.
*/

template <typename simplifier>
std::pair<double, double> calibrate(std::string filename, size_t N, simplifier s, std::pair<double, double> interval, std::string tag = "Current")
{

	int points = 0;
	int lines = 0;
	double eps = 1, avg_points = -1;
	// first, exponentiall grow the second interval
	std::cout << "Calibration started for " << tag << std::endl;
	while (true)
	{
		std::vector<size_t> npts{};

		double m = interval.second;
		scan(filename, "garbage.tmp",
				 [&npts, m, &s](linestring ls, emit_function_type emit, bool first) {
					 npts.push_back(s(ls, m).size());
				 });
		assert(npts.size() != 0);
		avg_points = static_cast<double>(std::accumulate(npts.begin(), npts.end(), 0)) / npts.size();
		//	     std::cout << "Calibrate to " << avg_points << std::endl;
		if (avg_points > N)
		{
			interval.second *= 2;
		}
		else
		{
			break;
		}
		std::cout << std::setprecision(5) << tag << ":" << avg_points << "Grow it between [" << interval.first << ";" << interval.second << "]";
	}
	double m;
	const double convergence_threshold = 1e-10;
	while (fabs(interval.second - interval.first) > convergence_threshold)
	{
		std::vector<size_t> npts;
		npts.clear();
		m = (interval.first + interval.second) / 2;
		scan(filename, "garbage.tmp",
				 [&npts, m, &s](linestring ls, emit_function_type emit, bool first) {
					 npts.push_back(s(ls, m).size());
				 });
		assert(npts.size() != 0);
		avg_points = static_cast<double>(std::accumulate(npts.begin(), npts.end(), 0)) / npts.size();
		//	     std::cout << "Calibrate to " << avg_points << std::endl;
		if (avg_points < N)
		{
			interval.second = m;
		}
		else
		{
			interval.first = m;
		}
#if OVERWRITE_PROGRESS
		std::cout << "\33[2K\r";
#endif
		std::cout << std::setprecision(5) << tag << ":\t" << avg_points << "\t Search Interval: [" << interval.first << ";" << interval.second << "]";

#if OVERWRITE_PROGRESS
		std::cout.flush();
#else
		std::cout << std::endl;
#endif
	}
#if OVERWRITE_PROGRESS
	std::cout << std::endl;
#endif

	return std::make_pair(m, avg_points);
}

/*
Algorithm 2.3: Run Calibration

runcalibration ()
   Reads paramters from JSON and tunes and evaluates all relevant algorithms writing
   aggregate results to a log file and to a JSON output. Therefore, parameters are
   estimated for all algorithms using parameters supplied from JSON. Then all data
   is transformed with found parameters and statistics on the Fr�chet distance as well as
   on the number of points are collected.
*/

void runcalibration()
{
	std::string filename = cfg.get<std::string>("input");
	int N = cfg.get<int>("N");
	std::cout << std::setprecision(32);

	// Calibrate Douglas-Peucker [dp.first ~ epsilon, dp.second ~ avg. number of points]
	auto dp = calibrate(
			filename, N, [](linestring in, double eps) {
				return douglas_peucker(in, eps);
			},
			{-1, 1}, "Douglas Peucker");
	// Calibrate Persistence with Beta Pruning
	auto beta = calibrate(
			filename, N,
			[](linestring in, double eps) {
				return persistence_beta(in, eps);
			},
			{0, 10000}, "Persistence Beta");
	// Calibrate MRS in two steps (beta for kappa*N, eps for N using this beta)
	double kappa = cfg.get<double>("kappa");
	auto mrs_beta = calibrate(
			filename, static_cast<int>(kappa * N),
			[](linestring in, double beta) {
				return persistence_beta(in, beta);
			},
			{0, 10000}, "MRS Kappa's Beta");

	auto mrs_eps = calibrate(
			filename, N,
			[&mrs_beta](linestring in, double eps) {
				return persistence_mrs(in, mrs_beta.first, eps, cfg.get<int>("mrs_rounds"));
			},
			{0, 1}, "Kappas MRS Epsilon");

	// Calibrate SDS in two steps (beta for kappa*N, eps for N using this beta)

	auto sds_beta = mrs_beta; // re-use result
	auto sds_eps = calibrate(
			filename, N,
			[&sds_beta](linestring in, double eps) {
				return persistence_sds(in, sds_beta.first, eps, cfg.get<int>("sds_rounds"));
			},
			{0, 1}, "Kappas Epsilon");

	// Phase 3:  Transform with complete statistics gathering and result writing
	std::map<std::string, std::vector<double>> errors;
	std::map<std::string, std::vector<double>> lengths;
	scan(filename, cfg.get<std::string>("output"),
			 [&](linestring ls, emit_function_type emit, bool first) {
				 if (first)
					 emit({"ls_orig", "ls_pers", "ls_dp", "ls_sds", "ls_mrs", "e_pers", "e_dp", "e_sds", "e_mrs",
								 "n_orig", "n_pers", "n_dp", "n_sds", "n_mrs"});
				 auto ls_pers = persistence_beta(ls, beta.first);
				 auto ls_dp = douglas_peucker(ls, dp.first);
				 auto ls_sds = persistence_sds(ls, sds_beta.first, sds_eps.first, cfg.get<int>("sds_rounds"));
				 auto ls_mrs = persistence_mrs(ls, mrs_beta.first, mrs_eps.first, cfg.get<int>("mrs_rounds"));
				 // errors

				 errors["pers"].push_back(frechet_distance(ls, ls_pers));
				 errors["dp"].push_back(frechet_distance(ls, ls_dp));
				 errors["sds"].push_back(frechet_distance(ls, ls_sds));
				 errors["mrs"].push_back(frechet_distance(ls, ls_mrs));

				 lengths["pers"].push_back(ls_pers.size());
				 lengths["dp"].push_back(ls_dp.size());
				 lengths["sds"].push_back(ls_sds.size());
				 lengths["mrs"].push_back(ls_mrs.size());

				 emit({
						 as_str(ls),
						 as_str(ls_pers),
						 as_str(ls_dp),
						 as_str(ls_sds),
						 as_str(ls_mrs),

						 std::to_string(errors["pers"].back()),
						 std::to_string(errors["dp"].back()),
						 std::to_string(errors["sds"].back()),
						 std::to_string(errors["mrs"].back()),
						 std::to_string(ls.size()),
						 std::to_string(ls_pers.size()),
						 std::to_string(ls_dp.size()),
						 std::to_string(ls_sds.size()),
						 std::to_string(ls_mrs.size()),

				 });
			 });

	// Compute average error statistics
	for (const auto &e : errors)
	{
		const std::string &key = e.first;
		const std::vector<double> &val = e.second;
		const std::vector<double> &_n = lengths[key];
		auto mean = std::accumulate(val.begin(), val.end(), 0.0) / val.size();
		auto meanLen = std::accumulate(_n.begin(), _n.end(), 0.0) / _n.size();
		std::cout << key << "\t:" << mean << "\t (" << meanLen << ")" << std::endl;
		cfg.set("e_" + key, mean);
		cfg.set("n_" + key, meanLen);
	}
}

std::string toString(std::vector<double> vec)
{
	std::stringstream ss;
	ss << "[";
	for (auto i = 0; i < vec.size(); i++)
	{
		if (i != 0)
			ss << ",";
		ss << vec[i];
	}
	ss << "]";
	return ss.str();
}
void profile(std::string filename, double pers_beta, double dp_eps, double sds_beta, double sds_eps, int sds_rounds, double mrs_beta, double mrs_eps, int mrs_rounds)
{
	int multiplier = 1;
	int iterations = 0;
	std::vector<persistence::Trajectory> trajectories;
	// extract all linestrings > 100 < 200
	scan(filename, cfg.get<std::string>("output"),
			 [&](linestring ls, emit_function_type emit, bool first) {
				 auto traj = extract(ls);
				 //printf("%d\n", traj.size());
				 trajectories.push_back(traj);
			 });

	std::sort(trajectories.begin(), trajectories.end(), [](const persistence::Trajectory &a, const persistence::Trajectory &b) { return a.size() < b.size(); });
	printf("size: %d min: %d max: %d\n", (int)trajectories.size(), (int)trajectories[0].size(), (int)trajectories.back().size());

	auto begin = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	double time;

	cfg.set("pers_time", (time / 1e+9) / (double)multiplier);

	std::vector<double> pers_times;
	for (auto j = 0; j < trajectories.size(); j++)
	{
		begin = std::chrono::high_resolution_clock::now();
		for (auto i = 0; i < multiplier; i++)
			persistence_beta(extract_linestring(trajectories[i]), pers_beta);

		end = std::chrono::high_resolution_clock::now();
		time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
		pers_times.push_back((time / 1e+9) / multiplier);
	}
	cfg.set("pers_times", toString(pers_times));
	printf("persistence_beta profiled!");
	std::vector<double> dp_times;
	for (auto j = 0; j < trajectories.size(); j++)
	{
		begin = std::chrono::high_resolution_clock::now();
		for (auto i = 0; i < multiplier; i++)
			douglas_peucker(extract_linestring(trajectories[i]), dp_eps);

		end = std::chrono::high_resolution_clock::now();
		time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
		dp_times.push_back((time / 1e+9) / multiplier);
	}

	cfg.set("dp_times", toString(dp_times));
	printf("douglas_peucker beta profiled!");

	std::vector<double> sds_times;
	for (auto j = 0; j < trajectories.size(); j++)
	{
		begin = std::chrono::high_resolution_clock::now();
		for (auto i = 0; i < multiplier; i++)
			persistence_sds(extract_linestring(trajectories[i]), sds_beta, sds_eps, sds_rounds);

		end = std::chrono::high_resolution_clock::now();
		time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
		sds_times.push_back((time / 1e+9) / multiplier);
	}

	cfg.set("sds_times", toString(sds_times));
	printf("persistence_sds profiled!");

	std::vector<double> mrs_times;
	for (auto j = 0; j < trajectories.size(); j++)
	{
		begin = std::chrono::high_resolution_clock::now();
		for (auto i = 0; i < multiplier; i++)
			persistence_mrs(extract_linestring(trajectories[i]), mrs_beta, mrs_eps, mrs_rounds);

		end = std::chrono::high_resolution_clock::now();
		time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
		mrs_times.push_back((time / 1e+9) / multiplier);
	}

	cfg.set("mrs_times", toString(mrs_times));
	printf("persistence_mrs profiled!");
}
void runprofile()
{
	std::string filename = cfg.get<std::string>("input");

	double pers_beta = cfg.get<double>("pers_beta");
	double dp_eps = cfg.get<double>("dp_eps");
	double sds_beta = cfg.get<double>("sds_beta");
	double sds_eps = cfg.get<double>("sds_eps");
	int sds_rounds = (int)cfg.get<double>("sds_rounds");
	double mrs_beta = cfg.get<double>("mrs_beta");
	double mrs_eps = cfg.get<double>("mrs_eps");
	int mrs_rounds = (int)cfg.get<double>("mrs_rounds");

	profile(filename, pers_beta, dp_eps, sds_beta, sds_eps, sds_rounds, mrs_beta, mrs_eps, mrs_rounds);
}

int main(int argc, char **argv)
{
	bool profile = false;
	// TODO retrieve profiling run from argc,argv
	if (argc < 3)
		throw(std::runtime_error("Call: " + std::string(argv[0]) + " <input> <output>, where input is a JSON and output as well. They are allowed to be in-place"));
	if (argc > 3)
		profile = std::string(argv[3]).compare("--profile") == 0;
	printf("%s", argv[3]);
	cfg.load(argv[1]);
	std::cout << cfg.dumps() << std::endl;
	if (!profile)
	{
		printf("run calibration: \n");
		runcalibration();
	}
	else
	{
		printf("run profile: \n");
		cfg.set("output", std::string(argv[2]));
		runprofile();
	}

	cfg.save(argv[2]);
	return 0;
}

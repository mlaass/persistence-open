#ifndef FRECHET_BOOST_GEOM_HPP
#define FRECHET_BOOST_GEOM_HPP

/*
  as well some stub methods to deal with std::vector<vec<double>> vs. linestring<point>
*/

std::vector<std::vector<double>> extract(linestring ls)
{
	std::vector<std::vector<double>> traj;
	std::transform(ls.begin(), ls.end(), std::back_inserter(traj),
								 [](point p) {
									 std::vector<double> ret{bg::get<0>(p), bg::get<1>(p)};
									 return ret;
								 });
	return traj;
}

linestring extract_linestring(const persistence::Result &result)
{
	linestring simplified;
	for (const auto &c : result.pruned)
	{
		simplified.push_back(point(c.vertex[0], c.vertex[1]));
	}
	return simplified;
}
linestring extract_linestring(std::vector<std::vector<double>> &result)
{
	linestring simplified;
	for (const auto &c : result)
	{
		simplified.push_back(point(c[0], c[1]));
	}
	return simplified;
}

linestring random_subset(linestring ls, size_t N)
{
	std::set<size_t> indices;
	if (N >= ls.size())
		throw(std::runtime_error("Cannot sample N from less than N"));
	while (indices.size() != N)
	{
		indices.insert(random() % ls.size());
	}
	linestring ret;
	for (const auto &i : indices)
		ret.push_back(ls[i]);
	return ret;
}

/*
provides adapters to use frechet distance with
boost geometry linestrings.

provides
double frechet_distance(const trajectory &t1, const trajectory &t2, int maxgrowth=1000)
*/

/*Frechet stub methods*/
struct get_coordinate
{
	template <size_t dim>
	static double get(const point &p)
	{
		return p.get<dim>();
	}
};

struct squared_distance
{
	double operator()(const point &p, const point &q) const
	{
		return bg::comparable_distance(p, q);
	}
};

template <typename trajectory>
double frechet_distance(const trajectory &t1, const trajectory &t2, int maxgrowth = 1000)
{
	static frechetrange::detail::dv::frechet_distance<2, get_coordinate,
																										squared_distance>
			fd;

	std::pair<double, double> interval = {0, 5};
	// grow right bound
	while (!fd.is_bounded_by(t1, t2, interval.second))
	{
		interval.first = interval.second;		// we know it is larger
		interval.second *= interval.second; // try this one
		if (--maxgrowth < 0)
			throw(std::runtime_error("Failed to grow Frechet distance"));
	}
	// now refine

	while ((interval.second - interval.first) > 0.001)
	{
		double m = (interval.first + interval.second) / 2;

		if (fd.is_bounded_by(t1, t2, m))
		{
			interval.second = m;
		}
		else
		{
			interval.first = m;
		}
	}
	/*	std::cout << "Final Interval: [" << interval.first << "," << interval.second
	<< "]" << std::endl;*/
	return ((interval.first + interval.second) / 2);
}

#endif

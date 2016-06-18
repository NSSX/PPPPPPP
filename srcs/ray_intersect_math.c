/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ray_intersect_math.c                               :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/03/15 10:30:55 by igomez            #+#    #+#             */
/*   Updated: 2015/03/15 17:01:26 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"
#include <math.h>

double		solve_2_degree(double a, double b, double c)
{
	double		tmp;
	double		res;

	if (fabs(a) < ERROR)
	{
		if (fabs(b) < ERROR)
			return (-1);
		res = -c / b;
	}
	else
	{
		tmp = b * b - 4 * a * c;
		if (tmp < 0)
			return (-1);
		tmp = sqrt(tmp);
		res = (-b - tmp) / (2 * a);
		if (res < 0)
			res = (-b + tmp) / (2 * a);
	}
	return (res);
}

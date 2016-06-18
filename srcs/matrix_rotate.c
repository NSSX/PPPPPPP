/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   matrix_rotate.c                                    :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/12/10 18:31:05 by igomez            #+#    #+#             */
/*   Updated: 2015/03/13 16:21:18 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1_matrix.h"
#include <math.h>

int			mat_rotate_x(t_matrix *buf, double ang, int sign)
{
	double		c;
	double		s;

	if (!buf || !sign)
		return (-1);
	c = cos(ang);
	s = sin((sign > 0) ? ang : -ang);
	*buf = M_NEW(4, 4);
	M_IJ(buf, 0, 0) = 1;
	M_IJ(buf, 1, 1) = c;
	M_IJ(buf, 1, 2) = -s;
	M_IJ(buf, 2, 1) = s;
	M_IJ(buf, 2, 2) = c;
	M_IJ(buf, 3, 3) = 1;
	return (0);
}

int			mat_rotate_y(t_matrix *buf, double ang, int sign)
{
	double		c;
	double		s;

	if (!buf || !sign)
		return (-1);
	c = cos(ang);
	s = sin((sign > 0) ? ang : -ang);
	*buf = M_NEW(4, 4);
	M_IJ(buf, 0, 0) = c;
	M_IJ(buf, 0, 2) = s;
	M_IJ(buf, 1, 1) = 1;
	M_IJ(buf, 2, 0) = -s;
	M_IJ(buf, 2, 2) = c;
	M_IJ(buf, 3, 3) = 1;
	return (0);
}

int			mat_rotate_z(t_matrix *buf, double ang, int sign)
{
	double		c;
	double		s;

	if (!buf || !sign)
		return (-1);
	c = cos(ang);
	s = sin((sign > 0) ? ang : -ang);
	*buf = M_NEW(4, 4);
	M_IJ(buf, 0, 0) = c;
	M_IJ(buf, 0, 1) = -s;
	M_IJ(buf, 1, 0) = s;
	M_IJ(buf, 1, 1) = c;
	M_IJ(buf, 2, 2) = 1;
	M_IJ(buf, 3, 3) = 1;
	return (0);
}

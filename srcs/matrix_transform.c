/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   matrix_transform.c                                 :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/12/10 18:31:05 by igomez            #+#    #+#             */
/*   Updated: 2015/03/15 12:52:20 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1_matrix.h"
#include <math.h>

int		mat_translate(t_matrix *buf, t_matrix *dir, int sign)
{
	if (dir->m != 4 || dir->n != 1 || !sign)
		return (-1);
	*buf = M_NEW(4, 4);
	if (sign > 0)
	{
		M_IJ(buf, 0, 3) = dir->mat[0];
		M_IJ(buf, 1, 3) = dir->mat[1];
		M_IJ(buf, 2, 3) = dir->mat[2];
	}
	else
	{
		M_IJ(buf, 0, 3) = -dir->mat[0];
		M_IJ(buf, 1, 3) = -dir->mat[1];
		M_IJ(buf, 2, 3) = -dir->mat[2];
	}
	M_IJ(buf, 0, 0) = 1;
	M_IJ(buf, 1, 1) = 1;
	M_IJ(buf, 2, 2) = 1;
	M_IJ(buf, 3, 3) = 1;
	return (0);
}

int		mat_scale(t_matrix *buf, t_matrix *dir, int sign)
{
	if (dir->m != 4 || dir->n != 1 || !sign)
		return (-1);
	*buf = M_NEW(4, 4);
	if (sign > 0)
	{
		M_IJ(buf, 0, 0) = dir->mat[0];
		M_IJ(buf, 1, 1) = dir->mat[1];
		M_IJ(buf, 2, 2) = dir->mat[2];
	}
	else
	{
		M_IJ(buf, 0, 0) = 1 / dir->mat[0];
		M_IJ(buf, 1, 1) = 1 / dir->mat[1];
		M_IJ(buf, 2, 2) = 1 / dir->mat[2];
	}
	M_IJ(buf, 3, 3) = 1;
	return (0);
}

int		mat_normalize(t_matrix *buf)
{
	double		norm;

	if (buf->m != 4 || buf->n != 1)
		return (-1);
	norm = sqrt(M_IJ(buf, 0, 0) * M_IJ(buf, 0, 0)
				+ M_IJ(buf, 1, 0) * M_IJ(buf, 1, 0)
				+ M_IJ(buf, 2, 0) * M_IJ(buf, 2, 0));
	M_IJ(buf, 0, 0) /= norm;
	M_IJ(buf, 1, 0) /= norm;
	M_IJ(buf, 2, 0) /= norm;
	return (0);
}

int		mat_setup(t_matrix *buf, t_matrix *pos, t_matrix *ang, t_matrix *scale)
{
	t_matrix	transfo;
	t_matrix	tmp;

	mat_translate(buf, pos, -1);
	mat_rotate_z(&transfo, M_IJ(ang, 2, 0), -1);
	mat_prod(&tmp, &transfo, buf);
	*buf = tmp;
	mat_rotate_y(&transfo, M_IJ(ang, 1, 0), -1);
	mat_prod(&tmp, &transfo, buf);
	*buf = tmp;
	mat_rotate_x(&transfo, M_IJ(ang, 0, 0), -1);
	mat_prod(&tmp, &transfo, buf);
	*buf = tmp;
	mat_scale(&transfo, scale, -1);
	mat_prod(&tmp, &transfo, buf);
	*buf = tmp;
	return (0);
}

int		mat_cleanup(t_matrix *buf, t_matrix *pos, t_matrix *ang)
{
	t_matrix	transfo;
	t_matrix	tmp;

	mat_rotate_x(buf, M_IJ(ang, 0, 0), 1);
	mat_rotate_y(&transfo, M_IJ(ang, 1, 0), 1);
	mat_prod(&tmp, &transfo, buf);
	*buf = tmp;
	mat_rotate_z(&transfo, M_IJ(ang, 2, 0), 1);
	mat_prod(&tmp, &transfo, buf);
	*buf = tmp;
	mat_translate(&transfo, pos, 1);
	mat_prod(&tmp, &transfo, buf);
	*buf = tmp;
	return (0);
}

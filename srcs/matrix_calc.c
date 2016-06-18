/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   matrix_calc.c                                      :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/12/02 13:13:13 by igomez            #+#    #+#             */
/*   Updated: 2016/06/16 17:54:10 by avella           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1_matrix.h"
#include "math.h"

int			mat_prod(t_matrix *buf, const t_matrix *m1, const t_matrix *m2)
{
	int		i;
	int		j;
	int		k;

	if (!buf || !m1 || !m2 || (m1->n != m2->m))
		return (-1);
	*buf = M_NEW(m1->m, m2->n);
	i = 0;
	while (i < buf->m)
	{
		j = 0;
		while (j < buf->n)
		{
			k = 0;
			while (k < m1->n)
			{
				M_IJ(buf, i, j) += (M_IJ(m1, i, k) * M_IJ(m2, k, j));
				k++;
			}
			j++;
		}
		i++;
	}
	return (0);
}

double		mat_dot_product(t_matrix *m1, t_matrix *m2)
{
	return (M_IJ(m1, 0, 0) * M_IJ(m2, 0, 0)
	+ M_IJ(m1, 1, 0) * M_IJ(m2, 1, 0)
	+ M_IJ(m1, 2, 0) * M_IJ(m2, 2, 0));
}

double		mat_norm(t_matrix *m1)
{
	return (sqrt((M_IJ(m1, 0, 0) * M_IJ(m1, 0, 0)
	+ M_IJ(m1, 1, 0) * M_IJ(m1, 1, 0)
	+ M_IJ(m1, 2, 0) * M_IJ(m1, 2, 0))));
}

double		mat_dist(t_matrix *m1, t_matrix *m2)
{
	t_matrix	tmp;

	tmp = M_DIR(
		M_IJ(m1, 0, 0) - M_IJ(m2, 0, 0),
		M_IJ(m1, 1, 0) - M_IJ(m2, 1, 0),
		M_IJ(m1, 2, 0) - M_IJ(m2, 2, 0));
	return (mat_norm(&tmp));
}

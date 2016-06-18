/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   color.c                                            :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/03/11 13:39:21 by igomez            #+#    #+#             */
/*   Updated: 2016/06/16 16:40:28 by avella           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"

int			color_scale(int c, double d)
{
	int		col[3];
	int		i;

	if (d < 0)
		return (c);
	col[0] = (c >> 16) & 0xff;
	col[1] = (c >> 8) & 0xff;
	col[2] = (c >> 0) & 0xff;
	i = 0;
	while (i < 3)
	{
		col[i] = d * col[i];
		if (col[i] > 255)
			col[i] = 255;
		++i;
	}
	return ((col[0] << 16) | (col[1] << 8) | (col[2]));
}

int			color_add(int a, int b)
{
	int		col[3];
	int		i;

	col[0] = ((a >> 16) & 0xff) + ((b >> 16) & 0xff);
	col[1] = ((a >> 8) & 0xff) + ((b >> 8) & 0xff);
	col[2] = ((a >> 0) & 0xff) + ((b >> 0) & 0xff);
	i = 0;
	while (i < 3)
	{
		if (col[i] > 255)
			col[i] = 255;
		++i;
	}
	return ((col[0] << 16) | (col[1] << 8) | (col[2]));
}

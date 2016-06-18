/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   debug.c                                            :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/03/05 15:38:08 by igomez            #+#    #+#             */
/*   Updated: 2015/03/18 13:49:43 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"
#include <stdio.h>

void		ft_printmat(t_matrix *buf)
{
	int		i;
	int		j;

	i = 0;
	while (i < buf->m)
	{
		j = 0;
		while (j < buf->n)
		{
			printf("%7.3f ", M_IJ(buf, i, j));
			j++;
		}
		printf("\n");
		i++;
	}
	printf("\n");
	fflush(stdout);
}

void		ft_printobj(t_obj *obj)
{
	printf("=== OBJ ===\n\n");
	switch (obj->param.type)
	{
		case CONE :
		printf("type:\nCONE\n\n");
		break ;
		case SPHERE :
		printf("type:\nSPHERE\n");
		break ;
		case CYLINDER :
		printf("type:\nCYLINDER\n\n");
		break ;
		case PLANE :
		printf("type:\nPLANE\n\n");
		break ;
	}
	printf("position:\n");
	ft_printmat(&obj->position);
	printf("ang:\n");
	ft_printmat(&obj->ang);
	printf("scale:\n");
	ft_printmat(&obj->scale);
	printf("color   : %#x\n", (unsigned int)obj->param.color);
	printf("setup:\n");
	ft_printmat(&obj->setup);
	printf("cleanup:\n");
	ft_printmat(&obj->cleanup);
	fflush(stdout);
}

void		ft_printray(t_ray *ray)
{
	printf("=== RAY ===\n\n");
	printf("start    :\n");
	ft_printmat(&ray->start);
	printf("direction:\n");
	ft_printmat(&ray->dir);
	if (ray->obj)
	{
		printf("intersection    :\n");
		ft_printmat(&(ray->intersection));
		printf("normal:\n");
		ft_printmat(&ray->normal);
	}
	fflush(stdout);
}

void		ft_printelapsed(t_env *e)
{
	printf("Calculation time: %f\n",
		(double)(e->c_new - e->c_old) / (double)CLOCKS_PER_SEC);
	fflush(stdout);
}

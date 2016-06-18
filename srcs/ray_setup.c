/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ray_setup.c                                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/03/07 13:15:44 by igomez            #+#    #+#             */
/*   Updated: 2016/06/17 18:37:05 by avella           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"

int			ray_setup_camera(t_ray *ray, t_cam *cam)
{
	t_matrix tmp;

	tmp = ray->start;
	mat_prod(&(ray->start), &(cam->setup), &tmp);
	tmp = ray->dir;
	mat_prod(&(ray->dir), &(cam->setup), &tmp);
	ray->obj = NULL;
	ray->dist = 0;
	return (0);
}

int			ray_setup_obj(t_ray *ray, t_obj *obj)
{
	t_matrix tmp;

	tmp = ray->start;
	mat_prod(&(ray->start), &(obj->setup), &tmp);
	tmp = ray->dir;
	mat_prod(&(ray->dir), &(obj->setup), &tmp);
	ray->obj = NULL;
	ray->dist = 0;
	return (0);
}

int			ray_cleanup_obj(t_ray *ray, t_obj *obj)
{
	t_matrix tmp;

	if (!ray->obj)
		return (-1);
	M_IJ(&(ray->intersection), 0, 0) *= M_IJ(&(obj->scale), 0, 0);
	M_IJ(&(ray->intersection), 1, 0) *= M_IJ(&(obj->scale), 1, 0);
	M_IJ(&(ray->intersection), 2, 0) *= M_IJ(&(obj->scale), 2, 0);
	tmp = ray->intersection;
	mat_prod(&(ray->intersection), &(obj->cleanup), &tmp);
	M_IJ(&(ray->normal), 0, 0) /= M_IJ(&(obj->scale), 0, 0);
	M_IJ(&(ray->normal), 1, 0) /= M_IJ(&(obj->scale), 1, 0);
	M_IJ(&(ray->normal), 2, 0) /= M_IJ(&(obj->scale), 2, 0);
	mat_normalize(&(ray->normal));
	tmp = ray->normal;
	mat_prod(&(ray->normal), &(obj->cleanup), &tmp);
	return (0);
}

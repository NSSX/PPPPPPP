/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   recalc_scene.c                                     :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/03/05 17:48:56 by igomez            #+#    #+#             */
/*   Updated: 2016/06/18 15:17:03 by avella           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct                                          s_vec3d
{
        double x;
        double y;
        double z;
}                                                       t_vec3d;

double          dot_product(t_vec3d *a, t_vec3d *b)
{
        return ((a->x * b->x) + (a->y * b->y) + (a->z * b->z));
}

void           normalize(t_vec3d *vec)
{
        double n;

        n = 1.0 / sqrt((vec->x * vec->x) + (vec->y * vec->y) + (vec->z * vec->z));
        vec->x *= n;
        vec->y *= n;
        vec->z *= n;
}

double          lim(double x, double min, double max)
{
        if (x < min)
                x = min;
        else if (x > max)
                x = max;
        return (x);
}

void            limit_for_vec(t_vec3d *vec, double a, double b)
{
        vec->x = lim(vec->x, a, b);
        vec->y = lim(vec->y, a, b);
        vec->z = lim(vec->z, a, b);
}

t_vec3d  ft_reflex(t_vec3d *incident, t_vec3d *n)
{
  t_vec3d  v;

  v.x = incident->x - 2.0 * dot_product(n, incident) * n->x;
  v.y = incident->y - 2.0 * dot_product(n, incident) * n->y;
  v.z = incident->z - 2.0 * dot_product(n, incident) * n->z;
  return (v);
}

int give_color(t_ray *ray, t_env *e)
{
	int color;

	ray_intersect(ray, e);
	color = ray_color(ray, e);
	if(e->nb_reflexion < 1 && ray->obj && ray->obj->param.type == PLANE)
	{

		e->nb_reflexion++;
		t_vec3d reflex;
		t_vec3d incident;
		t_vec3d normal;
		t_vec3d pos;

            pos.x =  M_IJ(&(ray->intersection), 0, 0);
            pos.y =  M_IJ(&(ray->intersection), 1, 0);
            pos.z =  M_IJ(&(ray->intersection), 2, 0);
			incident = (t_vec3d){M_IJ(&(ray->dir), 0, 0),M_IJ(&(ray->dir), 1, 0),M_IJ(&(ray->dir), 2, 0)};
            normal = (t_vec3d){M_IJ(&(ray->normal), 0, 0),
				M_IJ(&(ray->normal), 1, 0),
			M_IJ(&(ray->normal), 2, 0)};


            reflex = ft_reflex(&incident, &normal);
			normalize(&reflex);

            /*normal.x *= 0.001;
			normal.y *= 0.001;
			normal.z *= 0.001;*/
			pos.x += normal.x;
			pos.y += normal.y;
			pos.z += normal.z;

            M_IJ(&(ray->start), 0, 0) = pos.x;
			M_IJ(&(ray->start), 1, 0) = pos.y;
			M_IJ(&(ray->start), 2, 0) = pos.z;
			M_IJ(&(ray->dir), 0, 0) = reflex.x;
			M_IJ(&(ray->dir), 1, 0) = reflex.y;
			M_IJ(&(ray->dir), 2, 0) = reflex.z;

            int		col[3];

        	col[0] = (color >> 16) & 0xff;
        	col[1] = (color >> 8) & 0xff;
        	col[2] = (color >> 0) & 0xff;

            col[0] *= 0.25;
            col[1] *= 0.25;
            col[2] *= 0.25;


            int color_tmp = give_color(ray,e);

            int col_tmp[3];

            col_tmp[0] = (color_tmp >> 16) & 0xff;
        	col_tmp[1] = (color_tmp >> 8) & 0xff;
        	col_tmp[2] = (color_tmp >> 0) & 0xff;

            col_tmp[0] *= 1 - 0.25;
            col_tmp[1] *= 1 - 0.25;
            col_tmp[2] *= 1 - 0.25;

            col[0] += col_tmp[0];
            col[1] += col_tmp[1];
            col[2] += col_tmp[2];
            color = (col[0] << 16) | (col[1] << 8) | (col[2]);
    	   }
	return (color);
}

int		recalc_scene(t_env *e)
{
	int		x;
	int		y;
	t_ray	ray;

	e->c_old = clock();
	y = -1;
	while (++y < WIN_HEIGHT)
	{
		x = -1;
		while (++x < WIN_WIDTH)
		{
			e->nb_reflexion = 0;
			ray.start = M_POINT(0, 0, 0);
			ray.dir = M_DIR(DIST_TO_PROJ,
							WIN_WIDTH / 2 - x,
							WIN_HEIGHT / 2 - y);
			mat_normalize(&(ray.dir));
			ray_setup_camera(&ray, &(e->scene->cam));
			ft_draw_pix(TPIX(x, y, give_color(&ray,e)), e->screen);
		}
	}
	e->c_new = clock();
	return (0);
}

/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   env_scene_init.c                                   :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/03/03 16:15:06 by igomez            #+#    #+#             */
/*   Updated: 2016/06/18 15:11:22 by avella           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"
#include <math.h>

static int		env_scene_init_cam(t_scene *s)
{
	s->cam = (t_cam){M_POINT(-35, 9, 2), M_DIR(0, 0, 0), M_NEW(4, 4)};
	mat_cleanup(&(s->cam.setup), &(s->cam.position), &(s->cam.ang));
	return (0);
}

static int		env_scene_init_obj_add(t_obj *obj, t_scene *s)
{
	mat_setup(&obj->setup, &obj->position, &obj->ang, &obj->scale);
	mat_cleanup(&obj->cleanup, &obj->position, &obj->ang);
	ft_lstadd(&s->obj, ft_lstnew(obj, sizeof(t_obj)));
	return (0);
}

static int		env_scene_init_obj(t_scene *s)
{
	t_obj		tmp_obj;

	tmp_obj.position = M_POINT(0, 0, -10);
	tmp_obj.ang = ANG_PLANE_Z;
	tmp_obj.scale = M_DIR(1, 1, 1);
	tmp_obj.param = PARAM_OBJ_DEFAULT(PLANE, 0x00B2EE);
	env_scene_init_obj_add(&tmp_obj, s);

	tmp_obj.position = M_POINT(9, 2.5, -1);
	tmp_obj.ang = M_DIR(M_PI / 2, 0, 0);
	tmp_obj.scale = M_DIR(1, 1, 1);
	tmp_obj.param = PARAM_OBJ_DEFAULT(SPHERE, 0xff0000);
	env_scene_init_obj_add(&tmp_obj, s);

	tmp_obj.position = M_POINT(8, -5, 2.5);
	tmp_obj.ang = M_DIR(0, 0, 0);
	tmp_obj.scale = M_DIR(2, 2, 2);
	tmp_obj.param = PARAM_OBJ_DEFAULT(SPHERE, 0x00ff00);
	env_scene_init_obj_add(&tmp_obj, s);

	tmp_obj.position = M_POINT(10, -5, 0);
	tmp_obj.ang = M_DIR(-M_PI / 6, M_PI, 0);
	tmp_obj.scale = M_DIR(1, 4, 1);
	tmp_obj.param = PARAM_OBJ_DEFAULT(CYLINDER, 0xFFA500);
	env_scene_init_obj_add(&tmp_obj, s);

	tmp_obj.position = M_POINT(50, 25, 10);
	tmp_obj.ang = M_DIR(M_PI / 6, 0, 0);
	tmp_obj.scale = M_DIR(3, 1, 3);
	tmp_obj.param = PARAM_OBJ_DEFAULT(CONE, 0xff0f39);
	env_scene_init_obj_add(&tmp_obj, s);
	return (0);
}

static int		env_scene_init_spot(t_scene *s)
{
	t_spot		tmp_spot;

	tmp_spot.position = M_POINT(-100, 50, 50);
	tmp_spot.param = PARAM_SPOT_DEFAULT;
	ft_lstadd(&s->spot, ft_lstnew(&tmp_spot, sizeof(t_spot)));
	tmp_spot.position = M_POINT(-100, -50, 50);
	tmp_spot.param = PARAM_SPOT_DEFAULT;
	ft_lstadd(&s->spot, ft_lstnew(&tmp_spot, sizeof(t_spot)));
	return (0);
}

int				env_scene_init(t_env *e)
{
	if (!(e->scene = new_scene()))
		return (-1);
	env_scene_init_cam(e->scene);
	env_scene_init_obj(e->scene);
	env_scene_init_spot(e->scene);
	return (0);
}

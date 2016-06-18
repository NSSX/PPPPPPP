/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   env_mlx_init.c                                     :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/12/31 16:38:09 by igomez            #+#    #+#             */
/*   Updated: 2015/03/18 14:36:29 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"
#include "mlx.h"
#include "X.h"
#include <stdlib.h>

t_image				*image_new(t_env *e, int width, int height)
{
	t_image		*res;

	if (!(res = (t_image *)ft_memalloc(sizeof(t_image))))
		return (NULL);
	ft_bzero(res, sizeof(t_image));
	if (!(res->img = mlx_new_image(e->mlx, width, height)))
		return (NULL);
	if (!(res->data = mlx_get_data_addr(res->img, &(res->bpp),
		&(res->size_line), &(res->endian))))
		return (NULL);
	res->width = width;
	res->height = height;
	res->bypp = res->bpp / 8;
	return (res);
}

static int			env_mlx_init_hook(t_env *e)
{
	mlx_hook(e->win, Expose, ExposureMask, expose_hook, e);
	mlx_hook(e->win, KeyPress, KeyPressMask, key_hook, e);
	mlx_loop_hook(e->mlx, loop_hook, e);
	return (0);
}

int					env_mlx_init(t_env *e)
{
	e->width = WIN_WIDTH;
	e->height = WIN_HEIGHT;
	if (!(e->title = ft_strdup(TITLE)))
		return (-1);
	if (!(e->mlx = mlx_init()))
		return (-1);
	if (!(e->win = mlx_new_window(e->mlx, e->width, e->height, e->title)))
		return (-1);
	if (!(e->screen = image_new(e, e->width, e->height)))
		return (-1);
	env_mlx_init_hook(e);
	return (0);
}

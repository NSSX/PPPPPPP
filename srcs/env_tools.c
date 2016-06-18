/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   env_tools.c                                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/12/31 16:44:47 by igomez            #+#    #+#             */
/*   Updated: 2015/03/18 13:14:04 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"
#include "mlx.h"
#include <stdlib.h>

void	clean_env(t_env *e)
{
	free(e->title);
	if (e->screen->img)
		mlx_destroy_image(e->mlx, e->screen->img);
	free(e->screen);
	if (e->win)
		mlx_destroy_window(e->mlx, e->win);
	free(e->mlx);
}

void	ft_exit(int i, t_env *e)
{
	clean_env(e);
	exit(i);
}

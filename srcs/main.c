/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   main.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/12/31 16:29:05 by igomez            #+#    #+#             */
/*   Updated: 2015/03/18 14:37:53 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rtv1.h"
#include "mlx.h"

int		main(void)
{
	t_env		e;

	ft_bzero(&e, sizeof(t_env));
	if (env_mlx_init(&e) < 0)
		ft_putendl_fd("RTv1: mlx: Unexpected error", 2);
	else if (env_scene_init(&e) < 0)
		ft_putendl_fd("RTv1: scene: Unexpected error", 2);
	else
	{
		e.flag = 0;
		FLAG_SET(e.flag, MASK_RECALC);
		mlx_loop(e.mlx);
	}
	return (0);
}

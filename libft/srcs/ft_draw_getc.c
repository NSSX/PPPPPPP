/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_draw_getc.c                                     :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/02/16 11:00:07 by igomez            #+#    #+#             */
/*   Updated: 2015/02/16 14:03:47 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

int			ft_draw_getc(int x, int y, int *buf, t_image *img)
{
	ft_memcpy(buf, img->data + img->bypp * (y * img->width + x), img->bypp);
	return (0);
}

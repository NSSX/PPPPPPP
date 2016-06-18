/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_draw_square.c                                   :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/02/10 12:12:37 by igomez            #+#    #+#             */
/*   Updated: 2015/02/17 14:20:44 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include <libft.h>

int		ft_draw_square(t_pix p, t_pix delta, t_image *img)
{
	int		i;
	int		j;

	i = 0;
	while (i < delta.y)
	{
		delta.c = ((p.y + i) * img->width + p.x) * img->bypp;
		j = 0;
		while (j < img->bypp * delta.x)
		{
			((unsigned char *)(img->data + delta.c))[j] =
				((unsigned char *)(&(p.c)))[j % img->bypp];
			j++;
		}
		i++;
	}
	return (0);
}

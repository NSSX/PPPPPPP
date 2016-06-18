/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_putnbr.c                                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/06 11:33:19 by igomez            #+#    #+#             */
/*   Updated: 2014/12/01 22:02:34 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

void	ft_putnbr(int i)
{
	int div;
	int temp;

	div = 1;
	if (i < 0)
		ft_putchar('-');
	while (!((i / div > -10) && (i / div < 10)))
		div *= 10;
	while (div > 0)
	{
		temp = (i < 0) ? -(i / div) : i / div;
		ft_putchar('0' + temp);
		i %= div;
		div /= 10;
	}
}

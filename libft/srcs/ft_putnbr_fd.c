/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_putnbr_fd.c                                     :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/06 11:49:03 by igomez            #+#    #+#             */
/*   Updated: 2014/12/01 22:02:48 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"

void	ft_putnbr_fd(int i, int fd)
{
	int div;
	int temp;

	div = 1;
	if (i < 0)
		ft_putchar_fd('-', fd);
	while (!((i / div > -10) && (i / div < 10)))
		div *= 10;
	while (div > 0)
	{
		temp = (i < 0) ? -(i / div) : i / div;
		ft_putchar_fd('0' + temp, fd);
		i %= div;
		div /= 10;
	}
}
